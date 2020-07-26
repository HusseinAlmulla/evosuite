package org.evosuite.ga.metaheuristics;

import java.util.Map;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.evosuite.utils.LoggingUtils;

/**
 * This is using RBF value function and softmax policy
 * 
 * @author Hussein
 *
 */
public class RBF_SoftMaxPolicy {

	private static double alpha = 0;
	private static double beta = 0;
	private static double averageReward = 0;

	static double[] varTotal = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };
	static double[] meanTotal = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };

	public Hashtable<Integer, double[]> Weights_Action = new Hashtable<Integer, double[]>();
	private static int numOfOption = 0;
	int[] actionList;

	public RBF_SoftMaxPolicy() {

		averageReward = 0.0;/// numOfTilings;
		beta = 0.1;
		alpha = 0.4;
		numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;

		actionList = new int[numOfOption];
		for (int i = 0; i < numOfOption; i++) {
			Weights_Action.put(i, new double[] { 0.1, 0.1, 0.1, 0.1, 0.1 });
			actionList[i] = i;
		}
	}

	/**
	 * calculate the variance based on Welford's online algorithm
	 * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
	 * 
	 * @param counter
	 * @param features
	 * @return
	 */
	public void CalculateNewVarianceAndMean(int counter, double[] features) {
		for (int i = 0; i < features.length; i++) {

			double newMeanTotal = (meanTotal[i] * (counter - 1) + features[i]) / counter;
			double step1 = (counter - 1) * varTotal[i];
			double step2 = (features[i] - newMeanTotal) * (features[i] - meanTotal[i]);
			varTotal[i] = (step1 + step2) / counter;
			meanTotal[i] = newMeanTotal;
		}
	}

	public double value(int counter, double coverage, double numOfCoveredGoal, double fitness, double FoundGoals,
			int action) {

		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals,
				(coverage + numOfCoveredGoal - fitness + FoundGoals) };

		if (counter < 10) {
			varTotal = new double[] { 1, 1, 1, 1, 1 };
			meanTotal = feature_s_a;
		} else
			CalculateNewVarianceAndMean(counter, feature_s_a);
		double[] weights = Weights_Action.get(action);
		double sum = 0;

		for (int i = 0; i < feature_s_a.length; i++) {
			// Calculate the RBF for each feature
			double step1 = Math.pow((feature_s_a[i] - meanTotal[i]), 2);
			double step2 = 2 * varTotal[i];
			double x_i = Math.exp(-1 * (step1 / step2));

			sum += x_i * weights[i];
		}
		return sum;
	}

	public int SoftMaxPolicy(int counter, double coverage, double numOfCoveredGoal, double fitness, double FoundGoals) {

		
//		int max = 0;
		double softMaxBeta = 0.15;
		double sumOfValueFucntionAllAction = 0;
		double[] proValueActions = new double[numOfOption];
//		double[] val1 = new double[numOfOption];
//		double[] exp1 = new double[numOfOption];
		
		for (int action = 0; action < numOfOption; action++) {
//			val1[action] = value(counter, coverage, numOfCoveredGoal, fitness, FoundGoals, action);
//			exp1[action] = Math.exp(softMaxBeta * val1[action]);
			sumOfValueFucntionAllAction += Math
					.exp(value(counter, coverage, numOfCoveredGoal, fitness, FoundGoals, action) / softMaxBeta);

//			if (String.valueOf(sumOfValueFucntionAllAction) == "Infinity") {
//				double[] w = Weights_Action.get(action);
////				LoggingUtils.getEvoLogger().info(" weight " + w[0] + " " + w[1] + " " + w[2] + " " + w[3] + " " + w[4] );
//				LoggingUtils.getEvoLogger().info(" Renormalizing ++++++++++++++++");
//				renormalize_weights();
////				LoggingUtils.getEvoLogger().info(" weight " + w[0] + " " + w[1] + " " + w[2] + " " + w[3] + " " + w[4] );
//				return action;
//			}
		}

//		LoggingUtils.getEvoLogger().info("sumOfValueFucntionAllAction " + sumOfValueFucntionAllAction);
		for (int action = 0; action < numOfOption; action++) {
			double val_a = Math
					.exp(value(counter, coverage, numOfCoveredGoal, fitness, FoundGoals, action) / softMaxBeta);
			proValueActions[action] = val_a / sumOfValueFucntionAllAction;
		}

		EnumeratedIntegerDistribution dist = new EnumeratedIntegerDistribution(actionList, proValueActions);
		return dist.sample();
	}

	/**
	 * The weight sometime become to big and resulting of infinity which lead to
	 * exception
	 */
	private void renormalize_weights() {

		for (int i = 0; i < numOfOption; i++) {
			double[] w = Weights_Action.get(i);
			for (int j = 0; j < w.length; j++) {
				w[j] = w[j] / 1000;
			}
			Weights_Action.put(i, w);
		}
	}

	public void learn(int counter, double coverage, double numOfCoveredGoal, double fitness, double FoundGoals,
			int action, double newcoverage, double newnumOfCoveredGoal, double newfitness, double newFoundGoals,
			int newaction, double reward) {

		double estimation = value(counter, coverage, numOfCoveredGoal, fitness, FoundGoals, action);
		double value = value(counter, newcoverage, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
		double delta = reward - averageReward + value - estimation;

		averageReward += beta * delta;
		delta *= alpha;

		double[] feature_s_a = new double[] { coverage, numOfCoveredGoal, fitness, FoundGoals,
				(coverage + numOfCoveredGoal - fitness + FoundGoals) };
		double[] res = vectorTimeInt(feature_s_a, delta);

		double[] weights = Weights_Action.get(action);
		for (int i = 0; i < weights.length; i++)
			weights[i] += res[i];

		Weights_Action.put(action, weights);

		// LoggingUtils.getEvoLogger().info(" VAR " + varTotal + " MEAN " + meanTotal);
	}

	private double[] vectorTimeInt(double[] features, double z) {
		double[] res = new double[features.length];
		for (int i = 0; i < features.length; i++) {
			res[i] = features[i] * z;
		}
		return res;
	}

}
