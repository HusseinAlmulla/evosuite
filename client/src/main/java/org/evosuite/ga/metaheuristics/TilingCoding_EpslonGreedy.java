package org.evosuite.ga.metaheuristics;

/*
Original source is in Python
https://github.com/ShangtongZhang/reinforcement-learning-an-introduction
*/
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

import org.evosuite.utils.LoggingUtils;

public class TilingCoding_EpslonGreedy {

	private static int numOfTilings = 0;
	private static double alpha = 0;
	private static double beta = 0;
	private static double averageReward = 0;
	IHT hashTable;
	private double[] weights;
	
	private static int numOfOption = 0;
	public TilingCoding_EpslonGreedy(int nOfT, double al, double be) {
		numOfTilings = nOfT;
		int maxSize = 4096;
		hashTable = new IHT(maxSize);		
		weights = new double [maxSize];		
        averageReward = 0.0;
        alpha = al ;/// numOfTilings;
        beta = be;		
        numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;
	}
	
	public ArrayList<Integer> getActiveTiles(double[] floats, int action) { // coverage, int size, int NumOfCoveredGoals,double fitness
		
		ArrayList<Integer> activeTiles = hashTable.tiles(numOfTilings, floats, action);
		return activeTiles;
	}
	
	public double value(double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals, int action) {
		double[] floats = {coverage, size, numOfCoveredGoal, fitness, newFoundGoals, (coverage + numOfCoveredGoal - fitness + newFoundGoals)/4000};
		//LoggingUtils.getEvoLogger().info("Coverage " + coverage + "   Size " +  size, "   Num  "+  numOfCoveredGoal, "  Fit  " +fitness);
		ArrayList<Integer> activeTiles = getActiveTiles(floats, action);
		double sum = 0;
		//String st = "";
		for(int i=0; i < activeTiles.size(); i++) {
			int[] features = hashTable.index_features.get(activeTiles.get(i));
			double res = vectorTimeInt(features, weights[activeTiles.get(i)]);
			sum += res;
		
			
			
//			sum += weights[activeTiles.get(i)];
		}	        
		//LoggingUtils.getEvoLogger().info("Tiles " + st);
	        return sum;
	}
	
	private double vectorTimeInt(int[] features, double weight) {
		double sum = 0;
		for (int v: features) {
			sum += v * weight;
		}
		return sum;
	}
	
	public double[] stateValue(double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals) {
		double[] values = new double[numOfOption];
		int max = 0;
		for(int i=0; i<numOfOption; i++) {
			int action = i;
			values[i] = value(coverage, size, numOfCoveredGoal, fitness, newFoundGoals, action);
			//LoggingUtils.getEvoLogger().info("i "  + i + " " + values[i]);
			
			if(values[i] > values[max])
				max = i;
		}
		//LoggingUtils.getEvoLogger().info("max "  + max + " " + values[max]);
		double[] val = {values[max],max};
		return val;
	}
	
	
//	public double stateValue_softmax(double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals, int current_action) {
//		double[] values = new double[numOfOption];
//		int max = 0;
//		double current_action_value = value(coverage, size, numOfCoveredGoal, fitness, newFoundGoals, current_action);
//		double sum = 0;
//		for(int i=0; i<numOfOption; i++) {
//			int action = i;
//			values[i] = value(coverage, size, numOfCoveredGoal, fitness, newFoundGoals, action);
//			sum += Math.exp(0.5 * values[i]);			
//		}
//		double pro_pi = Math.exp(current_action_value)/sum;		
//		return pro_pi;
//	}
//	
	
	
	private static ArrayList<Double> semi_gediant = new ArrayList<Double>();
	
	public void learn(double coverage, int size, int numOfCoveredGoal, double fitness, int FoundGoals, int action
					, double newcoverage, int newsize, int newnumOfCoveredGoal, double newfitness, int newFoundGoals, int newaction, double reward) {
		
		double[] floats = {coverage, size, numOfCoveredGoal, fitness,FoundGoals, (coverage + numOfCoveredGoal - fitness + newFoundGoals)};
		ArrayList<Integer> activeTiles = getActiveTiles(floats, action);
		
        double estimation = 0;
		for(int i=0; i < activeTiles.size(); i++) {
			int[] features = hashTable.index_features.get(activeTiles.get(i));
			double res = vectorTimeInt(features, weights[activeTiles.get(i)]);
			estimation += res;
			
//			estimation += weights[activeTiles.get(i)];	
		}	
		
		// for maximizing the reward
		double value = value(newcoverage, newsize, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
        double delta = reward - averageReward + value - estimation;
//        LoggingUtils.getEvoLogger().info("estimated  Tail value"  + estimation + "  Tail Value  " + value);
		//for minimizing the reward
		//double delta =   reward + averageReward + value(newcoverage, newsize, newnumOfCoveredGoal, newfitness, newaction) - estimation;
		
        averageReward += beta * delta;
        delta *= alpha;
        
        double res = 0;
        for (int activeTile: activeTiles) {
        	int[] features = hashTable.index_features.get(activeTile);
			res += vectorTimeInt(features,delta);
			
//        	weights[activeTile] += delta;
			weights[activeTile] += res;
        }
	}
	
	public void printHash() {
		Hashtable<int[], Integer> d = hashTable.dictionary;
		Set<int[]> keys = d.keySet();
		for(int[] key: keys){
			String st = "";
			for(int k:key)
				st += k +" ";
			LoggingUtils.getEvoLogger().info("key " + st + "value " + hashTable.dictionary.get(key));
		}
	}
		
		public double[] printStatesValues(double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals) {
			double[] values = new double[numOfOption];
			int max = 0;
			for(int i=0; i<numOfOption; i++) {
				int action = i;
				values[i] = value(coverage, size, numOfCoveredGoal, fitness, newFoundGoals, action);
				LoggingUtils.getEvoLogger().info("i "  + i + " " + values[i]);
				if(values[i] > values[max])
					max = i;
			}
			System.out.println("max " + max + " " + values[max]);
			double[] val = {values[max],max};
			return val;
	}
	
	
	
	
	
	
	
}



