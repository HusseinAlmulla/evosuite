/**
 * Copyright (C) 2010-2018 Gordon Fraser, Andrea Arcuri and EvoSuite
 * contributors
 *
 * This file is part of EvoSuite.
 *
 * EvoSuite is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3.0 of the License, or
 * (at your option) any later version.
 *
 * EvoSuite is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with EvoSuite. If not, see <http://www.gnu.org/licenses/>.
 */
package org.evosuite.coverage.diversity;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.evosuite.Properties;
import org.evosuite.Properties.Criterion;
import org.evosuite.coverage.method.MethodCoverageSuiteFitness;
import org.evosuite.testcase.ExecutableChromosome;
import org.evosuite.testcase.TestChromosome;
import org.evosuite.testcase.execution.ExecutionResult;
import org.evosuite.testcase.statements.ArrayStatement;
import org.evosuite.testcase.statements.EntityWithParametersStatement;
import org.evosuite.testcase.statements.FieldStatement;
import org.evosuite.testcase.statements.MethodStatement;
import org.evosuite.testcase.statements.NullStatement;
import org.evosuite.testcase.statements.PrimitiveStatement;
import org.evosuite.testcase.statements.Statement;
import org.evosuite.testcase.statements.StringPrimitiveStatement;
import org.evosuite.testcase.statements.numeric.*;
import org.evosuite.testcase.statements.numeric.NumericalPrimitiveStatement;
import org.evosuite.testcase.variable.VariableReference;
import org.evosuite.testsuite.AbstractTestSuiteChromosome;
import org.evosuite.testsuite.TestSuiteChromosome;
import org.evosuite.testsuite.TestSuiteFitnessFunction;
import org.evosuite.utils.LoggingUtils;
import org.objectweb.asm.Type;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import org.evosuite.testcase.*;
/**
 * Fitness function for a whole test suite for all methods, including exceptional behaviour
 *
 * @author Gordon Fraser, Jose Miguel Rojas
 */
public class DiversityCoverageSuiteFitness extends TestSuiteFitnessFunction {

	private static final long serialVersionUID = 7980724720746662786L;


	private final static Logger logger = LoggerFactory.getLogger(DiversityCoverageSuiteFitness.class);
	private static TestSuiteFitnessFunction methodFitnessFunction = null;

	private double diversityCoverage;

	public DiversityCoverageSuiteFitness() {
		//methodFitnessFunction = new MethodCoverageSuiteFitness();		
	}


	/**
	 * {@inheritDoc}
	 *
	 * Execute all tests and count covered branches
	 */
	



	
	
	public double getFitness(AbstractTestSuiteChromosome<? extends ExecutableChromosome> suite) {

		//list of testcase along side with their map of parameters
		List<HashMap<Integer,String>> listTestCase = new ArrayList<>();
		List<HashMap<String,String>> listTestCaseRev = new ArrayList<>();
		
		// map of the each statement with their parameter
		HashMap<Integer,String> stmtParam = new HashMap<>();
		HashMap<String,String> param = new HashMap<>();

		TestSuiteChromosome testSuiteCopy = (TestSuiteChromosome) suite.clone();
		for (TestCase testCase : testSuiteCopy.getTests()) {	
			param = new HashMap<>();

			for ( int indx=0; indx < testCase.size(); indx++) {
				// problem here we are saving all variable not just the parameters
				Statement stmt = testCase.getStatement(indx);
				String superClass = stmt.getClass().getGenericSuperclass().toString();

				if(superClass.contains("NumericalPrimitiveStatement")) {
					stmtParam.put(indx, ((NumericalPrimitiveStatement<?>) stmt).getValue().toString());
					param.put("int",((NumericalPrimitiveStatement<?>) stmt).getValue().toString());					
				}
				else if(stmt instanceof NullStatement) {
					stmtParam.put(indx, stmt.getCode());
					param.put("str",stmt.getCode().toString());
				}
				else if(stmt instanceof FieldStatement ) {
					stmtParam.put(indx, ((FieldStatement) stmt).getField().toString());
					param.put("str",((FieldStatement) stmt).getField().toString());
				}

				else if(stmt instanceof StringPrimitiveStatement) {
					stmtParam.put(indx, ((StringPrimitiveStatement) stmt).getValue().toString());
					param.put("str",((StringPrimitiveStatement) stmt).getValue().toString());
				}
				else if (stmt instanceof ArrayStatement) {
					stmtParam.put(indx, ((ArrayStatement) stmt).getArrayReference().toString());
					param.put("str",((ArrayStatement) stmt).getArrayReference().toString());
				}
				else if (superClass.contains("EntityWithParametersStatement")) {					 
					EntityWithParametersStatement st1 = (EntityWithParametersStatement) stmt;
					List<VariableReference> list = st1.getParameterReferences();
					String st = "";

					if(list.size() > 0) {
						for(VariableReference para : list) {
							String p = stmtParam.get(para.getStPosition());
							if (p == null) 
								continue;
							st = st.concat(p).concat(" ");							 
						}
						st = st.trim();
						stmtParam.put(indx, st1.getMethodName()+":"+st);
						param.put("str",st1.getMethodName().concat(":").concat(st));
					}
					else {
						stmtParam.put(indx, st1.getMethodName());
						param.put("str",st1.getMethodName());
					}
				}
				else {
					stmtParam.put(indx, String.valueOf(stmt.getCode().toString()));
					param.put("str",String.valueOf(stmt.getCode().toString()));
				}
			}						
			listTestCaseRev.add(param);
			listTestCase.add(stmtParam);
		}
		
		double fitness1 = calculateDiversity5(listTestCaseRev);
//		fitness1 = fitness1 == 0 ? Double.MIN_VALUE : fitness1;
//		Criterion[] tmp = Properties.CRITERION;
//		Properties.CRITERION = new Criterion[] { Criterion.DIVERSITY, Criterion.METHOD };	
//		MethodCoverageSuiteFitness methodSuiteFitness = new MethodCoverageSuiteFitness();
//		double fitness2 = methodSuiteFitness.getFitness(testSuiteCopy);
//		double fitness = -1.0 * (fitness1 + fitness2);
//		double fitness = -1.0 * fitness1;
		
		double fitness = 1.0 / (1.0 + fitness1);
		fitness = Math.round(fitness * 10000000.0) / 10000000.0;
		updateIndividual(this, suite, fitness);
		suite.setDistance(fitness1);
		suite.setCoverage(this, (diversityCoverage));
//		Properties.CRITERION = tmp;	

		return fitness1;
	}

	
	private double calculateDiversity5(List<HashMap<String,String>> listPram) {

		double distance = 0;
		int listSize = listPram.size();
		double counter = 0;
		if (listPram.size() <= 1 )	 return 0;
		else {
			listPram.add(listPram.get(0));	
			for (int i = 0 ; i < listSize; i++) {
				for (int j = i+1 ; j < listSize; j++) {
				
					HashMap<String,String> testCase1 = listPram.get(i);
					HashMap<String,String> testCase2 = listPram.get((j));

				counter++;
				for (Entry<String, String> test1 : testCase1.entrySet()) 
					for (Entry<String, String> test2 : testCase2.entrySet()) {
						double tmpDistance = LevenshteinStringDistance(test1.getValue(), test2.getValue());
//						double tmpDistance = LSH_HammingStringDistance(test1.getValue(), test2.getValue());		
						distance += tmpDistance;
					}
				}
			}
		}
		return distance;
	}


	private double LevenshteinStringDistance(String str1, String str2)
	{	
		int[][] dis = new int[str2.length()+1][str1.length()+1];

		for (int i = 0 ; i <= str2.length(); i++)
			dis[i][0] = i;

		for (int j = 0 ; j <= str1.length(); j++) 				
			dis[0][j] = j;

		for (int i = 1 ; i <= str2.length(); i++) {	
			for (int j = 1 ; j <= str1.length(); j++) {				
				int cost = str1.charAt(j-1) == str2.charAt(i-1) ? 0 : 1;			
				dis[i][j] = Math.min(dis[i-1][j] + 1 , Math.min(dis[i][j-1] + 1 , dis[i-1][j-1] + cost));
			}		
		}
		double strLen = (double)(str1.length() + str2.length());
		double res = dis[str2.length()][str1.length()];
		return res;
	}
	
	//Locality-Sensitive Hashing (LSH) with Hamming distance
	private double LSH_HammingStringDistance(String str1, String str2)
	{	
		
		String str1Bin = str1;// org.apache.commons.lang3.StringUtils.leftPad(Integer.toBinaryString(str1.hashCode()), 32, '0');
		String str2Bin = str2;// org.apache.commons.lang3.StringUtils.leftPad(Integer.toBinaryString(str2.hashCode()), 32, '0');
		
		if(str1Bin.length() > str2Bin.length())
			str2Bin = org.apache.commons.lang3.StringUtils.leftPad(str2Bin, str1Bin.length(), '0');
		if(str2Bin.length() > str1Bin.length())
			str1Bin = org.apache.commons.lang3.StringUtils.leftPad(str1Bin, str2Bin.length(), '0');
		
		int diff = 0;
		int totalSum = 0;
		int str1Sum = 0;
		int str2Sum = 0;
		
		for( int i=0 ; i < str1Bin.length(); i++) {
			totalSum +=  str1Bin.charAt(i) * str2Bin.charAt(i);
			str1Sum += Math.pow(str1Bin.charAt(i), 2);
			str2Sum += Math.pow(str2Bin.charAt(i), 2);
		}
		double res = (double)totalSum/(double)(Math.sqrt(str1Sum) * Math.sqrt(str2Sum));
		
		
		return Math.round(res * 100.0)/100.0;
	}

}
