package org.evosuite.ga.metaheuristics;

import java.util.Map;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.evosuite.utils.LoggingUtils;


/**
 * This is for the GSD-Sarsa 
 * @author Hussein
 *
 */
public class LinearCombination_SoftMaxPolicy {

	private static double alpha = 0;
	private static double beta = 0;
	private static double averageReward = 0;
	private static int index = 0;

	
	public Hashtable<Integer, int[]>  Features_StateAction_dec = new Hashtable<Integer, int[]>();
	public Hashtable<Integer, double[]>  Weights_Action = new Hashtable<Integer, double[]>();
	private static int numOfOption = 0;
	
	double[] weights = new double[] {0.1,0.1,0.1,0.1,0.1,0.1 };
	
	int[] actionList;

	public LinearCombination_SoftMaxPolicy(double al, double be) {

		
		averageReward = 0.0;
		beta = 0.3;//be;
		alpha = 0.4;//al;
		numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;

		for(int i=0; i < numOfOption; i++)
			Weights_Action.put(i, new double[] {0.01,0.01,0.01,0.1,0.1});
		 actionList = new int[numOfOption];
		 for (int action = 0; action < numOfOption; action++) {
			 actionList[action] = action;
		 }
	}

	
	public double value(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals, int action) {

		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals), action/100 };

		double sum = 0;
//		double[] weights = Weights_Action.get(action);
		
		for (int i = 0; i < feature_s_a.length; i++) {
			sum += feature_s_a[i] * weights[i];
		}
		return sum;
	}

	public int SoftMaxPolicy(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		double softMaxBeta = 0.3; 
		double sumOfValueFucntionAllAction = 0;
		double[] proValueActions = new double[numOfOption];
		double[] val1 = new double[numOfOption];
		double[] exp1 = new double[numOfOption];
		
 		for (int action = 0; action < numOfOption; action++) {	
			val1[action] = value(coverage, numOfCoveredGoal, fitness, FoundGoals, action);
			exp1[action] = Math.exp(softMaxBeta * val1[action]);
			sumOfValueFucntionAllAction += Math.exp(value(coverage, numOfCoveredGoal, fitness, FoundGoals, action) / softMaxBeta );
			if(Double.isInfinite(sumOfValueFucntionAllAction))
				LoggingUtils.getEvoLogger().info("proValueActions  " + proValueActions[action]);
		}
		
		for (int action = 0; action < numOfOption; action++) {
			double val_a = Math.exp(value(coverage, numOfCoveredGoal, fitness, FoundGoals, action) / softMaxBeta );
			proValueActions[action] = val_a / sumOfValueFucntionAllAction;
			
			if(Double.isNaN(proValueActions[action]))
				LoggingUtils.getEvoLogger().info("proValueActions  " + proValueActions[action]);
		}
				
		EnumeratedIntegerDistribution dist   = new EnumeratedIntegerDistribution(actionList, proValueActions);
		return dist.sample();
	}
	
	public int stateValue(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		double[] values = new double[numOfOption];
		int max = 0;
		for (int i = 0; i < numOfOption; i++) {
			int action = i;
			values[i] = value(coverage, numOfCoveredGoal, fitness, FoundGoals, action);
			if (values[i] > values[max])
				max = i;
		}
		return max;
	}
	
	
	public void learn(double coverage, double numOfCoveredGoal, double fitness, double FoundGoals, int action, double newcoverage,
			double newnumOfCoveredGoal, double newfitness, double newFoundGoals, int newaction, double reward) {
	
		double estimation = value(coverage, numOfCoveredGoal, fitness, FoundGoals, action);
		double value = value(newcoverage, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
		double delta = reward - averageReward + value - estimation;

		averageReward += beta * delta;
		delta *= alpha;
		
		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals), action/100 };
		double[] res = vectorTimeInt(feature_s_a, delta);
		
//		double[] weights = Weights_Action.get(action);
		for (int i = 0; i < weights.length; i++)
			weights[i] += res[i];

//		Weights_Action.put(action, weights);
	}

	private double[] vectorTimeInt(double[] features, double z) {
		double[] res = new double[features.length];
		for (int i=0; i < features.length; i++) {
			res[i]= features[i] * z;
		}
		return res;
	}

}







