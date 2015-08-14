/**
 * Copyright (C) 2011,2012 Gordon Fraser, Andrea Arcuri and EvoSuite
 * contributors
 * 
 * This file is part of EvoSuite.
 * 
 * EvoSuite is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 * 
 * EvoSuite is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Public License for more details.
 * 
 * You should have received a copy of the GNU Public License along with
 * EvoSuite. If not, see <http://www.gnu.org/licenses/>.
 */
package org.evosuite.idnaming;

import com.examples.with.different.packagename.idnaming.SimpleIdNaming;
import com.examples.with.different.packagename.sette.L4_Collections;
import com.examples.with.different.packagename.sette.SnippetInputContainer;
import com.examples.with.different.packagename.stable.Base64;
import com.examples.with.different.packagename.staticfield.UnstableAssertion;

import org.evosuite.EvoSuite;
import org.evosuite.Properties;
import org.evosuite.SystemTest;
import org.evosuite.ga.metaheuristics.GeneticAlgorithm;
import org.evosuite.strategy.TestGenerationStrategy;
import org.evosuite.testsuite.TestSuiteChromosome;
import org.junit.Assert;
import org.junit.Test;

public class TestIDNaming extends SystemTest {

	@Test
	public void testIDNamingOn() {
		EvoSuite evosuite = new EvoSuite();

		String targetClass = SnippetInputContainer.class.getCanonicalName();

		Properties.TARGET_CLASS = targetClass;
		Properties.ID_NAMING = true;
		Properties.JUNIT_TESTS = true;
		StringBuilder analysisCriteria = new StringBuilder();
        analysisCriteria.append(Properties.Criterion.METHOD); analysisCriteria.append(",");
        analysisCriteria.append(Properties.Criterion.OUTPUT); analysisCriteria.append(",");
        analysisCriteria.append(Properties.Criterion.BRANCH);
        Properties.ANALYSIS_CRITERIA = analysisCriteria.toString();
        
        Properties.CRITERION = new Properties.Criterion[3];
        Properties.CRITERION[0] = Properties.Criterion.METHOD;
        Properties.CRITERION[1] = Properties.Criterion.OUTPUT;
        Properties.CRITERION[2] = Properties.Criterion.BRANCH;

		String[] command = new String[] { "-generateSuite", "-class", targetClass };

		Object result = evosuite.parseCommandLine(command);
		GeneticAlgorithm<?> ga = getGAFromResult(result);
		TestSuiteChromosome best = (TestSuiteChromosome) ga.getBestIndividual();
		System.out.println("EvolvedTestSuite:\n" + best);

		int goals = TestGenerationStrategy.getFitnessFactories().get(0).getCoverageGoals().size(); // assuming single fitness function
		Assert.assertEquals("Wrong number of goals: ", 15, goals);
		Assert.assertEquals("Non-optimal coverage: ", 1d, best.getCoverage(), 0.001);
	}

	@Test
	public void testIDNamingOff() {
		EvoSuite evosuite = new EvoSuite();

		String targetClass = SimpleIdNaming.class.getCanonicalName();

		Properties.TARGET_CLASS = targetClass;
		Properties.ID_NAMING = false;
		Properties.JUNIT_TESTS = true;

		String[] command = new String[] { "-generateSuite", "-class", targetClass };

		Object result = evosuite.parseCommandLine(command);
		GeneticAlgorithm<?> ga = getGAFromResult(result);
		TestSuiteChromosome best = (TestSuiteChromosome) ga.getBestIndividual();
		System.out.println("EvolvedTestSuite:\n" + best);

		int goals = TestGenerationStrategy.getFitnessFactories().get(0).getCoverageGoals().size(); // assuming single fitness function
		Assert.assertEquals("Wrong number of goals: ", 7, goals);
		Assert.assertEquals("Non-optimal coverage: ", 1d, best.getCoverage(), 0.001);
	}

}
