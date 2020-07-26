package org.evosuite.ga.metaheuristics;

import java.util.Map;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.evosuite.utils.LoggingUtils;

/**
 * This is for the Semi_Gradient TD(0)
 * @author Hussein
 *
 */
public class LinearCombination_SoftMaxPolicy2 {

	private static double alpha = 0;
	

	public Hashtable<Integer, int[]>  Features_StateAction_dec = new Hashtable<Integer, int[]>();
	public Hashtable<Integer, double[]>  Weights_Action = new Hashtable<Integer, double[]>();
	
	double[] weights = new double[5];
	private static int numOfOption = 0;
	
	int[] actionList;

	public LinearCombination_SoftMaxPolicy2(double al, double be) {

		alpha = al;
		numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;

		for(int i=0; i < numOfOption; i++)
			Weights_Action.put(i, new double[] {0.1,0.1,0.1,0.1,0.1});
		 actionList = new int[numOfOption];
		 for (int action = 0; action < numOfOption; action++) {
			 actionList[action] = action;
		 }
	}

	
	public double value(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals) };
		
		double sum = 0;		
		for (int i = 0; i < feature_s_a.length; i++) {
			sum += feature_s_a[i] * weights[i];
		}
		return sum;
	}

	public int SoftMaxPolicy(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		double softMaxBeta = 0.2; 
		double sumOfValueFucntionAllAction = 0;
		double[] proValueActions = new double[numOfOption];
		double[] val1 = new double[numOfOption];
		double[] exp1 = new double[numOfOption];
		int max = 0;
		
 		for (int action = 0; action < numOfOption; action++) {	
			val1[action] = value(coverage, numOfCoveredGoal, fitness, FoundGoals);
			exp1[action] = Math.exp(softMaxBeta * val1[action]);
			sumOfValueFucntionAllAction += Math.exp(value(coverage, numOfCoveredGoal, fitness, FoundGoals) / softMaxBeta );
		}
		
		for (int action = 0; action < numOfOption; action++) {
			double val_a = Math.exp(value(coverage, numOfCoveredGoal, fitness, FoundGoals) / softMaxBeta );
			proValueActions[action] = val_a / sumOfValueFucntionAllAction;
			
			if (proValueActions[action] > proValueActions[max])
				max = action;
		}
		
		EnumeratedIntegerDistribution dist   = new EnumeratedIntegerDistribution(actionList, proValueActions);
		max = dist.sample();

		return max;
	}
	
	public double stateValue(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		double values = value(coverage, numOfCoveredGoal, fitness, FoundGoals);
		return values;
	}
	
	/**
	 * Learning based on the algorithm of Semi-Gradient TD(0)
	 * @param coverage
	 * @param numOfCoveredGoal
	 * @param fitness
	 * @param FoundGoals
	 * @param action
	 * @param newcoverage
	 * @param newnumOfCoveredGoal
	 * @param newfitness
	 * @param newFoundGoals
	 * @param reward
	 */
	
	public void learn(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals, double newcoverage,
			double newnumOfCoveredGoal, double newfitness, double newFoundGoals, double reward) {
	
		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals) };
		
		double value = value(coverage, numOfCoveredGoal, fitness, FoundGoals);
		double estimation = value(newcoverage, newnumOfCoveredGoal, newfitness, newFoundGoals);
		
		double diff = alpha * (reward + (0.3 * estimation) - value);
		double[] error = vectorTimeInt(feature_s_a, diff);
		
		for (int i = 0; i < weights.length; i++)
			weights[i] += error[i];
	}

	private double[] vectorTimeInt(double[] features, double z) {
		double[] res = new double[features.length];
		for (int i=0; i < features.length; i++) {
			res[i]= features[i] * z;
		}
		return res;
	}

}







