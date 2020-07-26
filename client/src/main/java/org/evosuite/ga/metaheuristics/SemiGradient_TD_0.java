package org.evosuite.ga.metaheuristics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import org.evosuite.utils.LoggingUtils;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class SemiGradient_TD_0 {

	private static int numberOfOption;
	int[] numbers_of_selections;
	

	public SemiGradient_TD_0(int num) {
		numberOfOption = num;
		numbers_of_selections = new int[numberOfOption];
		init();

		/*
		 * try { wrt = new PrintWriter(new FileOutputStream(new File("D:\\trace.txt"),
		 * true)); } catch (FileNotFoundException e1) { e1.printStackTrace(); }
		 */
	}

	final List<Integer> shuffled_options = new ArrayList<Integer>();

	private void init() {

		for (int j = 0; j < numberOfOption; j++) {
			shuffled_options.add(j);
		}
		Collections.shuffle(shuffled_options);
		}


	// ValueFunction2 value_function;
	LinearCombination_SoftMaxPolicy2 value_function = new LinearCombination_SoftMaxPolicy2(0.5, 0.1);

	private double current_coverage;
	private int current_numOfCoveredGoal;
	private double current_fitness;
	private int current_foundGoal;
//	ArrayList<Integer> op = new ArrayList<Integer>();

	public int getCurrent_foundGoal() {
		return current_foundGoal;
	}

	public void setCurrent_foundGoal(int current_foundGoal) {
		this.current_foundGoal = current_foundGoal;
	}

	Random rand = new Random();

	
	private int getAction_softmax(int counter, double coverage, double numOfCoveredGoal, double fitness,
			double newFoundGoals) {

		int action = -1;
		if (counter < numberOfOption) {
			action = shuffled_options.get(counter);
			numbers_of_selections[action] += 1;
		} else {
//			if (counter % step == 0) {
//				step = (step > 2) ? (step - 1) : 2;
//				action = value_function.stateValue(coverage , numOfCoveredGoal, fitness, newFoundGoals);
//				selected[action] += 1;
//			} else {
				action = value_function.SoftMaxPolicy(coverage , numOfCoveredGoal, fitness, newFoundGoals);
				numbers_of_selections[action] += 1;
//			}
		}
		return action;
	}

	public int getCurrent_action(int counter) {
		
		double nc = (double) Math.round(current_coverage * 1000) / 1000;
		double nf = (double) Math.round(current_fitness) / 1000;
		double ncg = (double) Math.round(current_numOfCoveredGoal) / 1000;
		double nfg = (double) Math.round(current_foundGoal) / 1000;
		
		int newAction = getAction_softmax(counter, nc, ncg, nf, nfg);
		numbers_of_selections[newAction] +=1;
		return newAction;
	}

	double[] rwd;

	public void DSGSarsa_part2(int counter, double newCoverage, int newSize, int newNumOfCoveredGoal, double newFitness,
			int newFoundGoals, double reward_score) {

//		This used with linear combination of feature 
		double cc = (double) Math.round(current_coverage * 1000) / 1000;
		double nc = (double) Math.round(newCoverage * 1000) / 1000;
		
		double cf = (double) Math.round(current_fitness) / 1000;
		double nf = (double) Math.round(newFitness) / 1000;
		
		double ccg = (double) Math.round(current_numOfCoveredGoal) / 1000;
		double ncg = (double) Math.round(newNumOfCoveredGoal) / 1000;
		
		double cfg = (double) Math.round(current_foundGoal) / 1000;
		double nfg = (double) Math.round(newFoundGoals) / 1000;
		
		value_function.learn( cc, ccg, cf, cfg, nc, ncg, nf, nfg, reward_score);
		
		current_coverage = newCoverage;
		current_numOfCoveredGoal = newNumOfCoveredGoal;
		current_foundGoal = newFoundGoals;
		current_fitness = newFitness;
	}

	// return the final estimation for the best option that the system recognize
	public int DSGSarsa_FinalEstimation(double newCoverage, int newSize, int newNumOfCoveredGoal, double newFitness,
			int newFoundGoals) {

		double nc = (double) Math.round(newCoverage * 1000) / 1000;
		double nf = (double) Math.round(newFitness * 100) / 100;
		double ncg = (double) Math.round(newNumOfCoveredGoal * 1000) / 1000;
		double nfg = (double) Math.round(newFoundGoals * 1000) / 1000;
		
		int val = value_function.SoftMaxPolicy(nc , ncg,  nf,	nfg);

		LoggingUtils.getEvoLogger().info("  max " + val);
		return  val;
	}
	

	public void prt() {
		for (int i = 0; i < numberOfOption; i++) {
			LoggingUtils.getEvoLogger()
					.info("selected[" + i + "] " + numbers_of_selections[i]);
		}
	}

}
