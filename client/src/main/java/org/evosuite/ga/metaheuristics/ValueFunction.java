package org.evosuite.ga.metaheuristics;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

import org.evosuite.utils.LoggingUtils;

public class ValueFunction {

	private static int numOfTilings = 0;
	private static double alpha = 0;
	private static double beta = 0;
	private static double averageReward = 0;
	IHT hashTable;
	private double[] weights;
	
	private static int numOfOption = 0;
	public ValueFunction(int nOfT, double al, double be) {
		numOfTilings = nOfT;
		int maxSize = 4096;
		hashTable = new IHT(maxSize);		
		weights = new double [maxSize];		
        averageReward = 0.0;
        alpha = al / numOfTilings;
        beta = be;		
        numOfOption = org.evosuite.ga.metaheuristics.MonotonicGA.numberOfOption;
	}
	
	public ArrayList<Integer> getActiveTiles(double[] floats, int action) { // coverage, int size, int NumOfCoveredGoals,double fitness
		
		ArrayList<Integer> activeTiles = hashTable.tiles(numOfTilings, floats, action);
		return activeTiles;
	}
	
	public double value(double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals, int action) {
		double[] floats = {coverage, size, numOfCoveredGoal, fitness, newFoundGoals};
		//LoggingUtils.getEvoLogger().info("Coverage " + coverage + "   Size " +  size, "   Num  "+  numOfCoveredGoal, "  Fit  " +fitness);
		ArrayList<Integer> activeTiles = getActiveTiles(floats, action);
		double sum = 0;
		//String st = "";
		for(int i=0; i < activeTiles.size(); i++) {
			sum += weights[activeTiles.get(i)];
			//st += activeTiles.get(i) + " , ";
		}	        
		//LoggingUtils.getEvoLogger().info("Tiles " + st);
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
	
	
	public void learn(double coverage, int size, int numOfCoveredGoal, double fitness, int FoundGoals, int action
					, double newcoverage, int newsize, int newnumOfCoveredGoal, double newfitness, int newFoundGoals, int newaction, double reward) {
		
		double[] floats = {coverage, size, numOfCoveredGoal, fitness};
		ArrayList<Integer> activeTiles = getActiveTiles(floats, action);
        double estimation = 0;
		for(int i=0; i < activeTiles.size(); i++) {
			estimation += weights[activeTiles.get(i)];
		}	
		// for maximizing the reward
		double value = value(newcoverage, newsize, newnumOfCoveredGoal, newfitness, newFoundGoals, newaction);
        double delta = reward - averageReward + value - estimation;
        LoggingUtils.getEvoLogger().info("estimated  Tail value"  + estimation + "  Tail Value  " + value);
		//for minimizing the reward
		//double delta =   reward + averageReward + value(newcoverage, newsize, newnumOfCoveredGoal, newfitness, newaction) - estimation;
		
        averageReward += beta * delta;
        delta *= alpha;
        
        for (int activeTile: activeTiles)
            weights[activeTile] += delta;
        
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



