package org.evosuite.ga.metaheuristics;


import java.util.ArrayList;
import java.util.Random;
import org.evosuite.utils.LoggingUtils;


public class RLAlgorithm {

   public RLAlgorithm() {
	   init();
	   /*try {
			wrt = new PrintWriter(new FileOutputStream(new File("D:\\trace.txt"), true));
		}
		catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}*/
   }
   

	

	//////////////////// DSGSarsa ////////////////////////////
	
	ValueFunction value_function;
	private void init() {
		/*Random rnd = new Random();
		int sum = 0;
		for(int i=1; i<8; i++) {
			w[0][i] = rnd.nextInt(2);
			sum += w[0][i];
		}		
		for(int i=1; i<8; i++)
			w[0][i] = w[0][i] / sum;*/
		
		value_function = new ValueFunction(20, 0.1, 0.1);
	}
	
	//private double[][] w  = new double[92][8]; //{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
	
	/*private double R_average = 0;
	private double beta = 0.312;
	private double alpha = 0.12;
	private double delta = 0;
	
	private double[] q = new double[92];*/
	
	private int current_action = 0;
	private double current_coverage;
	private int current_size;
	private int current_numOfCoveredGoal;
	private double current_fitness;
	private int current_foundGoal;
	ArrayList<Integer> op = new ArrayList<Integer>();

	public int getCurrent_foundGoal() {
		return current_foundGoal;
	}

	public void setCurrent_foundGoal(int current_foundGoal) {
		this.current_foundGoal = current_foundGoal;
	}

	private int init = 1;
	
	Random rand = new Random();
	private int[] selected = new int[92]; 
	
	private int getAction(int counter, double coverage, int size, int numOfCoveredGoal, double fitness, int newFoundGoals) {
		
		int action = -1;		
		if(counter<92) {
			action = counter;
			selected[action] += 1;
			value_function.stateValue(coverage, size, numOfCoveredGoal, fitness, newFoundGoals);	
		}
		else {
			if(counter % 2 == 0) {
				double[] values = value_function.stateValue(coverage, size, numOfCoveredGoal, fitness, newFoundGoals);				
				action = (int) values[1];
				selected[action] += 1; 
			}
			else {
				action = Math.abs(rand.nextInt(92));
				selected[action] += 1;
			}
		}
		op.add(action);
		return action;
}
	
	public int getCurrent_action() {
		return current_action;
	}

	public void setCurrent_action(int current_action) {
		this.current_action = current_action;
	}
	double[] rwd = new double[92];
	
	public void DSGSarsa_part2(int counter, double newCoverage, int newSize, int newNumOfCoveredGoal, double newFitness,int newFoundGoals, double reward_score) {
		
		int newAction = getAction(counter, newCoverage, newSize, newNumOfCoveredGoal, newFitness, newFoundGoals);
		value_function.learn(current_coverage, current_size, current_numOfCoveredGoal, current_fitness, current_foundGoal, current_action, 
				newCoverage, newSize, newNumOfCoveredGoal, newFitness,newFoundGoals, newAction, reward_score);
		rwd[getCurrent_action()] += reward_score; 
		current_coverage = newCoverage;
		current_size = newSize;
		current_numOfCoveredGoal = newNumOfCoveredGoal; 
		current_fitness = newFitness;
		setCurrent_action(newAction);
	}	
    	

// return the final estimation for the best option that the system recognize 
	public int DSGSarsa_FinalEstimation(double newCoverage, int newSize, int newNumOfCoveredGoal, double newFitness, int newFoundGoals) {
		
		double[] val =value_function.stateValue(newCoverage, newSize, newNumOfCoveredGoal, newFitness, newFoundGoals);
		//double[] val =value_function.printStatesValues(newCoverage, newSize, newNumOfCoveredGoal, newFitness, newFoundGoals);
		LoggingUtils.getEvoLogger().info("state value  " + val[0] + "  max " + val[1]);
		return (int)val[1];		
		
	}
	
	/*public int DSGSarsa_part1(int counter) {
		int ad = -1;

		if(init<92) {
			ad = init;
			init++;
		}
		else {
			if(counter % 2 == 0) {
				ad = getMax();
				 selected[ad] += 1; 
			}
			else {
				ad = Math.abs(rand.nextInt(92));
				wrt.println("Random " + ad);
				selected[ad] += 1;
			}
				
		}
		
		wrt.println("iter "+ iter +"\nad- " + ad);
		Criterion[] cc = getCriteria(ad);
		String ss = "";
		for(int i=0; i<cc.length; i++)
			ss += cc[i].toString();
		wrt.println("critera  " + ss);
		
		cs = ns;
		wcs = wns;
		
		previous_action = current_action;
		current_action = ad;
		
		return ad;
	}*/
	
	/*public void DSGSarsa_part2(double reward_score) {
		delta = reward_score - R_average + q[current_action] - q[previous_action];
		R_average += beta * delta;
		
		wrt.println("Delta "+ delta +",  R_average " + R_average);
		String str = "";
		
		for(int i=0; i<8; i++) {
			w[current_action][i] = w[previous_action][i] +  (alpha * delta);
			str += w[i] +",";
		}
		wrt.println("W [" + str + "]");
	}	*/
	/*
	private int getMax() {
		int max = 0;
		for(int i=0; i< 92; i++) {
			double tmp = 0;
			for(int j=0; j<8; j++) {
				 tmp += w[i][j];// * x_s_a[i][j];
			}
			q[i] = tmp;
			if(q[i] > q[max])
				max = i;
		}
		return max;
	}*/
	
	/*private static Criterion[] mainCriteria = new Criterion[] {
            Criterion.LINE, Criterion.BRANCH, Criterion.EXCEPTION, Criterion.WEAKMUTATION, Criterion.OUTPUT, Criterion.METHOD, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH  };

	public Criterion[] getCriteria(int ad) {
		
		int[] SA = x_s_a[ad];
		int sum = IntStream.of(SA).sum();
		Criterion[] subCriteria = new Criterion[sum];
		int j =0;
		for(int i=0; i< SA.length; i++) {
			if(SA[i] == 1) {
				subCriteria[j] = mainCriteria[i];
				j++;
			}
		}
		return subCriteria;
	}
	
	*/
	
	public void prt() {
		int ms = 0;
		int mr = 0;
		for(int i=0; i< 92; i++) {
			//for(int j=0; j< 8; j++)
				LoggingUtils.getEvoLogger().info("selected[" + i + "] " + selected[i] + " -  reward  " + rwd[i]);
				if(selected[i] > selected[ms])
					ms = i;
				if(rwd[i] > rwd[mr])
					mr = i;							
		}
		//value_function.printHash();
		LoggingUtils.getEvoLogger().info(" Max selected " + ms + " -  max reward  " + mr);
		/*try {
			PrintWriter wrt = new PrintWriter(new FileOutputStream(new File("D:\\trace.txt"), true));
			for(int i=0; i<op.size(); i++)
				wrt.write(op.get(i) + " - ");
			wrt.close();			
		}
		catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}*/
		
	}

		
	
}
