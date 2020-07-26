package org.evosuite.ga.metaheuristics;

import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import org.evosuite.utils.LoggingUtils;

public class LinearCombination_EpslonGreedy {

	private static double alpha = 0;
	private static double beta = 0;
	private static double averageReward = 0;
	private static int index = 0;

//	Map<int[], int[]> Features_StateAction = new HashMap<int[], int[]>();
	public Hashtable<Integer, int[]>  Features_StateAction_dec = new Hashtable<Integer, int[]>();
	public Hashtable<Integer, double[]>  Weights_Action = new Hashtable<Integer, double[]>();
	private static int numOfOption = 0;
	

	public LinearCombination_EpslonGreedy(double al, double be) {

		averageReward = 0.0;
		alpha = al;/// numOfTilings;
		beta = be;
		numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;

		for(int i=0; i < numOfOption; i++)
			Weights_Action.put(i, new double[] {0.1, 0.1,0.1,0.1,0.1,0.1});
		
	}

	//////////////// NOT USED //////////////////////
	public int[] getFeatures(int states_action) { 

		int[] feature_s_a = (int[]) Features_StateAction_dec.get(states_action);
		return feature_s_a;
	}
	public void addFeatures_state_action( int coverage,int sizq, int numOfCoveredGoal, int fitness, int FoundGoals, int action) { 
		
		Integer states_action = (String.valueOf(coverage) + String.valueOf(numOfCoveredGoal) + String.valueOf(fitness) + String.valueOf(FoundGoals) + String.valueOf(action)).hashCode();
		int[] feature = { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals) };
		
		Features_StateAction_dec.put(states_action, feature);
	}
	/////////////////////////////////////
	
	
	public double value(int coverage, int size, int numOfCoveredGoal, int fitness, int FoundGoals, int action) {

//		Integer states_action = coverage + numOfCoveredGoal + fitness + FoundGoals + action;
		
//		Integer states_action = (String.valueOf(coverage)+ String.valueOf(size) + String.valueOf(numOfCoveredGoal) + String.valueOf(fitness) + String.valueOf(FoundGoals) + String.valueOf(action)).hashCode();
//		int[] feature_s_a = getFeatures(states_action);
//		if (feature_s_a == null)
//			return Double.MIN_VALUE;
//		
		
		
		int[] feature_s_a = new int[] { coverage, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals) };
		
		double sum = 0;
		double[] weights = Weights_Action.get(action);
		for (int i = 0; i < feature_s_a.length; i++) {
			sum += feature_s_a[i] * weights[i];
		}
		return sum;
	}

	public double[] stateValue(int coverage, int size, int numOfCoveredGoal, int fitness, int FoundGoals) {

		double[] values = new double[numOfOption];
		int max = 0;
		for (int i = 0; i < numOfOption; i++) {
			int action = i;
			values[i] = value(coverage, size, numOfCoveredGoal, fitness, FoundGoals, action);
			if (values[i] > values[max])
				max = i;
		}
		double[] val = { values[max], max };
		return val;
	}

	
	public double learn(int coverage, int size, int numOfCoveredGoal, int fitness, int FoundGoals, int action, 
			int newcoverage, int newSize, int newnumOfCoveredGoal, int newfitness, int newFoundGoals, int newaction, double reward) {

//		Integer states_action = coverage + numOfCoveredGoal + fitness + FoundGoals + action;
//		Integer states_action = (String.valueOf(coverage) + String.valueOf(size) + String.valueOf(numOfCoveredGoal) + String.valueOf(fitness) + String.valueOf(FoundGoals) + String.valueOf(action)).hashCode();
//		int[] feature_s_a = getFeatures(states_action);
//		if (feature_s_a == null) {
//			addFeatures_state_action(coverage, size, numOfCoveredGoal, fitness, FoundGoals, action);
//			feature_s_a = getFeatures(states_action);
//		}
		
		int[] feature_s_a = new int[] { coverage, size, numOfCoveredGoal, fitness, FoundGoals, (coverage + numOfCoveredGoal - fitness + FoundGoals) };
		
		double estimation = 0;
		double[] weights = Weights_Action.get(action);
		for (int i = 0; i < feature_s_a.length; i++) {
			estimation += feature_s_a[i] * weights[i];
		}

		double value = value(newcoverage, newSize, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
//		if (value == Double.MIN_VALUE) {
//			addFeatures_state_action(newcoverage, newSize, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
//			value = value(newcoverage, newSize, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
//		}
		double delta = reward - averageReward + value - estimation;

		averageReward += beta * delta;
		delta *= alpha;
				
		double[] res = vectorTimeInt(feature_s_a, delta);
		for (int i = 0; i < weights.length; i++)
			weights[i] += res[i];
		
		Weights_Action.put(action, weights);
		
		return value;		
	}

	private double[] vectorTimeInt(int[] features, double z) {
		double[] res = new double[features.length];
		for (int i=0; i < features.length; i++) {
			res[i]= features[i] * z;
		}
		return res;
	}

}
