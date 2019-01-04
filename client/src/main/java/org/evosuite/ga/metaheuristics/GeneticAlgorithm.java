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
package org.evosuite.ga.metaheuristics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.evosuite.Properties;
import org.evosuite.Properties.Algorithm;
import org.evosuite.Properties.Criterion;
import org.evosuite.coverage.CoverageCriteriaAnalyzer;
import org.evosuite.coverage.FitnessFunctions;
import org.evosuite.coverage.TestFitnessFactory;
import org.evosuite.ga.Chromosome;
import org.evosuite.ga.ChromosomeFactory;
import org.evosuite.ga.FitnessFunction;
import org.evosuite.ga.archive.Archive;
import org.evosuite.ga.bloatcontrol.BloatControlFunction;
import org.evosuite.ga.localsearch.DefaultLocalSearchObjective;
import org.evosuite.ga.localsearch.LocalSearchBudget;
import org.evosuite.ga.localsearch.LocalSearchObjective;
import org.evosuite.ga.operators.crossover.CrossOverFunction;
import org.evosuite.ga.operators.crossover.SinglePointCrossOver;
import org.evosuite.ga.operators.ranking.RankBasedPreferenceSorting;
import org.evosuite.ga.operators.ranking.RankingFunction;
import org.evosuite.ga.operators.selection.RankSelection;
import org.evosuite.ga.operators.selection.SelectionFunction;
import org.evosuite.ga.populationlimit.IndividualPopulationLimit;
import org.evosuite.ga.populationlimit.PopulationLimit;
import org.evosuite.ga.stoppingconditions.MaxGenerationStoppingCondition;
import org.evosuite.ga.stoppingconditions.StoppingCondition;
import org.evosuite.runtime.Random;
import org.evosuite.symbolic.DSEStats;
import org.evosuite.testcase.TestFitnessFunction;
import org.evosuite.testcase.execution.ExecutionTracer;
import org.evosuite.testsuite.TestSuiteChromosome;
import org.evosuite.testsuite.TestSuiteFitnessFunction;
import org.evosuite.utils.ArrayUtil;
import org.evosuite.utils.LoggingUtils;
import org.evosuite.utils.Randomness;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Abstract superclass of genetic algorithms
 * 
 * @author Gordon Fraser
 */
public abstract class GeneticAlgorithm<T extends Chromosome> implements SearchAlgorithm, Serializable {

	private static final long serialVersionUID = 5155609385855093435L;

	private static final Logger logger = LoggerFactory.getLogger(GeneticAlgorithm.class);

	/** Fitness function to rank individuals */
	protected List<FitnessFunction<T>> fitnessFunctions = new ArrayList<FitnessFunction<T>>();

	/** Selection function to select parents */
	protected SelectionFunction<T> selectionFunction = new RankSelection<T>();

	/** CrossOver function */
	protected CrossOverFunction crossoverFunction = new SinglePointCrossOver();

	/** Current population */
	protected List<T> population = new ArrayList<T>();

	/** Generator for initial population */
	protected ChromosomeFactory<T> chromosomeFactory;

	/** Listeners */
	protected transient Set<SearchListener> listeners = new HashSet<SearchListener>();

	/** List of conditions on which to end the search */
	protected transient Set<StoppingCondition> stoppingConditions = new HashSet<StoppingCondition>();

	/** Bloat control, to avoid too long chromosomes */
	protected Set<BloatControlFunction> bloatControl = new HashSet<BloatControlFunction>();

	/** Local search might need a different local objective */
	protected LocalSearchObjective<T> localObjective = new DefaultLocalSearchObjective<>();

	/** The population limit decides when an iteration is done */
	protected PopulationLimit populationLimit = new IndividualPopulationLimit();

	/** Age of the population */
	protected int currentIteration = 0;

	protected double localSearchProbability = Properties.LOCAL_SEARCH_PROBABILITY;

	/** Selected ranking strategy **/
	protected RankingFunction<T> rankingFunction = new RankBasedPreferenceSorting<T>();

	/**
	 * Constructor
	 * 
	 * @param factory
	 *            a {@link org.evosuite.ga.ChromosomeFactory} object.
	 */
	public GeneticAlgorithm(ChromosomeFactory<T> factory) {
		chromosomeFactory = factory;
		addStoppingCondition(new MaxGenerationStoppingCondition());
		if (Properties.LOCAL_SEARCH_RATE > 0)
			addListener(LocalSearchBudget.getInstance());
		// addBloatControl(new MaxSizeBloatControl());
	}

	/**
	 * Generate one new generation
	 */
	protected abstract void evolve();

	/**
	 * Local search is only applied every X generations
	 * 
	 * @return a boolean.
	 */
	protected boolean shouldApplyLocalSearch() {
		// If local search is not set to a rate, then we don't use it at all
		if (Properties.LOCAL_SEARCH_RATE <= 0)
			return false;

		if (getAge() % Properties.LOCAL_SEARCH_RATE == 0) {
			if (Randomness.nextDouble() <= localSearchProbability) {
				return true;
			}
		}
		return false;
	}

	protected void disableFirstSecondaryCriterion() {
		if (TestSuiteChromosome.getSecondaryObjectivesSize() > 1) {
			TestSuiteChromosome.disableFirstSecondaryObjective();
			if (ArrayUtil.contains(Properties.SECONDARY_OBJECTIVE, Properties.SecondaryObjective.IBRANCH)) {
				ExecutionTracer.disableContext();
			}
			logger.info("second secondary criterion enabled");
		}
	}

	protected void enableFirstSecondaryCriterion() {
		if (TestSuiteChromosome.getSecondaryObjectivesSize() > 1) {
			TestSuiteChromosome.enableFirstSecondaryObjective();
			if (ArrayUtil.contains(Properties.SECONDARY_OBJECTIVE, Properties.SecondaryObjective.IBRANCH)) {
				ExecutionTracer.enableContext();
			}
			logger.info("first secondary criterion enabled");
		}
	}

	/**
	 * enable and disable secondary criteria according to the strategy defined in
	 * the Properties file.
	 * 
	 * @param starvationCounter
	 */
	protected void updateSecondaryCriterion(int starvationCounter) {

		if (Properties.ENABLE_SECONDARY_OBJECTIVE_AFTER > 0 && TestSuiteChromosome.getSecondaryObjectivesSize() > 1) {

			double progress = this.progress() * 100.0;

			if (progress > Properties.ENABLE_SECONDARY_OBJECTIVE_AFTER) {
				if (Properties.ENABLE_SECONDARY_OBJECTIVE_STARVATION) {
					updateSecondaryObjectiveStarvation(starvationCounter);
				} else {
					enableFirstSecondaryCriterion();
					Properties.ENABLE_SECONDARY_OBJECTIVE_AFTER = 0;
				}
			}
		} else if (Properties.ENABLE_SECONDARY_OBJECTIVE_STARVATION && Properties.ENABLE_SECONDARY_OBJECTIVE_AFTER == 0
				&& TestSuiteChromosome.getSecondaryObjectivesSize() > 1) {
			updateSecondaryObjectiveStarvation(starvationCounter);
		}
	}

	private void updateSecondaryObjectiveStarvation(int starvationCounter) {
		if (starvationCounter > Properties.STARVATION_AFTER_GENERATION
				&& !TestSuiteChromosome.isFirstSecondaryObjectiveEnabled()) {
			enableFirstSecondaryCriterion();
		} else {
			if (starvationCounter == 0 && TestSuiteChromosome.isFirstSecondaryObjectiveEnabled()
					&& TestSuiteChromosome.getSecondaryObjectivesSize() > 1) {
				disableFirstSecondaryCriterion();
			}
		}
	}

	/**
	 * Apply local search, starting from the best individual and continue applying
	 * it to all individuals until the local search budget is used up.
	 * 
	 * The population list is re-ordered if needed.
	 */
	protected void applyLocalSearch() {
		if (!shouldApplyLocalSearch())
			return;

		logger.debug("Applying local search");
		LocalSearchBudget.getInstance().localSearchStarted();

		boolean improvement = false;

		for (Chromosome individual : population) {
			if (isFinished())
				break;

			if (LocalSearchBudget.getInstance().isFinished()) {
				logger.debug("Local search budget used up, exiting local search");
				break;
			}

			if (individual.localSearch(localObjective)) {
				improvement = true;
			}
		}

		if (improvement) {
			DSEStats.getInstance().reportNewIncrease();
			updateProbability(true);
			logger.debug("Increasing probability of applying LS to " + localSearchProbability);
		} else {
			DSEStats.getInstance().reportNewDecrease();
			updateProbability(false);
			logger.debug("Decreasing probability of applying LS to " + localSearchProbability);
		}

		if (improvement) {
			// If an improvement occurred to one of the individuals, it could
			// be the case that the improvement was so good, that the individual
			// has surpassed to the previous individual, which makes the population
			// list not sorted any more.
			if (!populationIsSorted()) {
				this.sortPopulation();
			}
		}
	}

	/**
	 * Returns true if the population is sorted according to the fitness values.
	 * 
	 * @return true if the population is sorted (or empty)
	 */
	private boolean populationIsSorted() {
		Chromosome previousIndividual = null;
		for (Chromosome currentIndividual : this.population) {
			if (previousIndividual != null) {
				if (!isBetterOrEqual(previousIndividual, currentIndividual)) {
					// at least two individuals are not sorted
					return false;
				}
			}
			previousIndividual = currentIndividual;
		}
		// the population is sorted (or empty)
		return true;
	}

	/**
	 * Returns true if the fitness functions are maximization functions or false if
	 * all fitness functions are minimisation functions. It expects that all
	 * fitnessFunctions are minimising or maximising, it cannot happen that
	 * minimization and maximization functions are together.
	 * 
	 * @return
	 */
	private boolean isMaximizationFunction() {
		return fitnessFunctions.get(0).isMaximizationFunction();
	}

	protected void updateProbability(boolean improvement) {
		if (improvement) {
			localSearchProbability *= Properties.LOCAL_SEARCH_ADAPTATION_RATE;
			localSearchProbability = Math.min(localSearchProbability, 1.0);
		} else {
			localSearchProbability /= Properties.LOCAL_SEARCH_ADAPTATION_RATE;
			localSearchProbability = Math.max(localSearchProbability, Double.MIN_VALUE);
		}
		// localSearchProbability = Math.pow(
		// 1.0 + ((1.0 - localSearchProbability) / localSearchProbability) *
		// Math.exp(delta), -1.0);
	}

	/**
	 * Apply dynamic symbolic execution
	 */
	/*
	 * @Deprecated protected void applyDSE() {
	 * logger.info("Applying DSE at generation " + currentIteration);
	 * DSEBudget.DSEStarted();
	 * 
	 * boolean success = false;
	 * 
	 * for (Chromosome individual : population) { if (isFinished()) break;
	 * 
	 * if (DSEBudget.isFinished()) break;
	 * 
	 * boolean result = individual.applyDSE(this); if(result) success = true; }
	 * 
	 * if(Properties.DSE_ADAPTIVE_PROBABILITY > 0.0) { if(success) {
	 * Properties.DSE_ADAPTIVE_PROBABILITY *= Properties.DSE_ADAPTIVE_RATE;
	 * Properties.DSE_ADAPTIVE_PROBABILITY =
	 * Math.min(Properties.DSE_ADAPTIVE_PROBABILITY, 1.0); } else {
	 * Properties.DSE_ADAPTIVE_PROBABILITY /= Properties.DSE_ADAPTIVE_RATE;
	 * Properties.DSE_ADAPTIVE_PROBABILITY =
	 * Math.max(Properties.DSE_ADAPTIVE_PROBABILITY, Double.MIN_VALUE); } logger.
	 * info("Updating DSE probability to "+Properties.DSE_ADAPTIVE_PROBABILITY); } }
	 */
	/**
	 * Set up initial population
	 */
	public abstract void initializePopulation();

	/**
	 * {@inheritDoc}
	 * 
	 * Generate solution
	 */
	@Override
	public abstract void generateSolution();

	/**
	 * Fills the population at first with recycled chromosomes - for more
	 * information see recycleChromosomes() and ChromosomeRecycler - and after that,
	 * the population is filled with random chromosomes.
	 * 
	 * This method guarantees at least a proportion of
	 * Properties.initially_enforeced_Randomness % of random chromosomes
	 * 
	 * @param population_size
	 *            a int.
	 */
	protected void generateInitialPopulation(int population_size) {
		generateRandomPopulation(population_size - population.size());
	}

	/**
	 * This method can be used to kick out chromosomes when the population is
	 * possibly overcrowded
	 * 
	 * Depending on the Property "starve_by_fitness" chromosome are either kicked
	 * out randomly or according to their fitness
	 * 
	 * @param limit
	 *            a int.
	 */
	protected void starveToLimit(int limit) {
		if (Properties.STARVE_BY_FITNESS)
			starveByFitness(limit);
		else
			starveRandomly(limit);
	}

	/**
	 * This method can be used to kick out random chromosomes in the current
	 * population until the given limit is reached again.
	 * 
	 * @param limit
	 *            a int.
	 */
	protected void starveRandomly(int limit) {
		while (population.size() > limit) {
			int removePos = Randomness.nextInt() % population.size();
			population.remove(removePos);
		}
	}

	/**
	 * This method can be used to kick out the worst chromosomes in the current
	 * population until the given limit is reached again.
	 * 
	 * @param limit
	 *            a int.
	 */
	protected void starveByFitness(int limit) {
		calculateFitnessAndSortPopulation();
		for (int i = population.size() - 1; i >= limit; i--) {
			population.remove(i);
		}
	}

	/**
	 * Generate random population of given size
	 * 
	 * @param population_size
	 *            a int.
	 */
	protected void generateRandomPopulation(int population_size) {
		logger.debug("Creating random population");
		for (int i = 0; i < population_size; i++) {
			T individual = chromosomeFactory.getChromosome();
			for (FitnessFunction<?> fitnessFunction : this.fitnessFunctions) {
				individual.addFitness(fitnessFunction);
			}

			population.add(individual);
			if (isFinished())
				break;
		}
		logger.debug("Created " + population.size() + " individuals");
	}

	/**
	 * Delete all current individuals
	 */
	public void clearPopulation() {
		logger.debug("Resetting population");
		population.clear();
	}

	/**
	 * Add new fitness function (i.e., for new mutation)
	 * 
	 * @param function
	 *            a {@link org.evosuite.ga.FitnessFunction} object.
	 */
	public void addFitnessFunction(FitnessFunction<T> function) {
		fitnessFunctions.add(function);
		localObjective.addFitnessFunction(function);
	}

	public void addFitnessFunctions(List<FitnessFunction<T>> functions) {
		for (FitnessFunction<T> function : functions)
			this.addFitnessFunction(function);
	}

	/**
	 * Get currently used fitness function
	 * 
	 * @return a {@link org.evosuite.ga.FitnessFunction} object.
	 */
	public FitnessFunction<T> getFitnessFunction() {
		return fitnessFunctions.get(0);
	}

	/**
	 * Get all used fitness function
	 * 
	 * @return a {@link org.evosuite.ga.FitnessFunction} object.
	 */
	public List<FitnessFunction<T>> getFitnessFunctions() {
		return fitnessFunctions;
	}

	public int getNumberOfFitnessFunctions() {
		return fitnessFunctions.size();
	}

	/**
	 * 
	 * @return
	 */
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();

		int i = 0;
		for (Chromosome c : population) {
			str.append("\n  - test " + i);

			for (FitnessFunction<T> ff : this.fitnessFunctions) {
				DecimalFormat df = new DecimalFormat("#.#####");
				str.append(", " + ff.getClass().getSimpleName().replace("CoverageSuiteFitness", "") + " "
						+ df.format(c.getFitness(ff)));
			}

			i++;
		}

		return str.toString();
	}

	/**
	 * Set new fitness function (i.e., for new mutation)
	 * 
	 * @param function
	 *            a {@link org.evosuite.ga.operators.selection.SelectionFunction}
	 *            object.
	 */
	public void setSelectionFunction(SelectionFunction<T> function) {
		selectionFunction = function;
	}

	/**
	 * Get currently used fitness function
	 * 
	 * @return a {@link org.evosuite.ga.operators.selection.SelectionFunction}
	 *         object.
	 */
	public SelectionFunction<T> getSelectionFunction() {
		return selectionFunction;
	}

	/**
	 * Set the new ranking function (only used by MOO algorithms)
	 * 
	 * @param function
	 *            a {@link org.evosuite.ga.operators.ranking.RankingFunction} object
	 */
	public void setRankingFunction(RankingFunction<T> function) {
		this.rankingFunction = function;
	}

	/**
	 * Get currently used ranking function (only used by MOO algorithms)
	 * 
	 * @return a {@link org.evosuite.ga.operators.ranking.RankingFunction} object
	 */
	public RankingFunction<T> getRankingFunction() {
		return this.rankingFunction;
	}

	/**
	 * Set new bloat control function
	 * 
	 * @param bloat_control
	 *            a {@link org.evosuite.ga.bloatcontrol.BloatControlFunction}
	 *            object.
	 */
	public void setBloatControl(BloatControlFunction bloat_control) {
		this.bloatControl.clear();
		addBloatControl(bloat_control);
	}

	/**
	 * Set new bloat control function
	 * 
	 * @param bloat_control
	 *            a {@link org.evosuite.ga.bloatcontrol.BloatControlFunction}
	 *            object.
	 */
	public void addBloatControl(BloatControlFunction bloat_control) {
		this.bloatControl.add(bloat_control);
	}

	/**
	 * Check whether individual is suitable according to bloat control functions
	 * 
	 * @param chromosome
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @return a boolean.
	 */
	public boolean isTooLong(Chromosome chromosome) {
		for (BloatControlFunction b : bloatControl) {
			if (b.isTooLong(chromosome))
				return true;
		}
		return false;
	}

	/**
	 * Get number of iterations
	 * 
	 * @return Number of iterations
	 */
	public int getAge() {
		return currentIteration;
	}

	/**
	 * Calculate fitness for all individuals
	 */
	protected void calculateFitness() {
		logger.debug("Calculating fitness for " + population.size() + " individuals");

		Iterator<T> iterator = this.population.iterator();
		while (iterator.hasNext()) {
			T c = iterator.next();
			if (isFinished()) {
				if (c.isChanged())
					iterator.remove();
			} else {
				this.calculateFitness(c);
			}
		}
	}

	/**
	 * Calculate fitness for an individual
	 * 
	 * @param c
	 */
	protected void calculateFitness(T c) {
		for (FitnessFunction<T> fitnessFunction : this.fitnessFunctions) {
			fitnessFunction.getFitness(c);
			notifyEvaluation(c);
		}
	}

	/**
	 * Calculate fitness for all individuals and sort them
	 */
	protected void calculateFitnessAndSortPopulation() {
		this.calculateFitness();
		// Sort population
		this.sortPopulation();
	}

	/**
	 * <p>
	 * getPopulationSize
	 * </p>
	 * 
	 * @return a int.
	 */
	public int getPopulationSize() {
		return population.size();
	}

	/**
	 * Copy best individuals
	 * 
	 * @return a {@link java.util.List} object.
	 */
	@SuppressWarnings("unchecked")
	protected List<T> elitism() {
		logger.debug("Elitism with ELITE = " + Properties.ELITE);

		List<T> elite = new ArrayList<T>();

		for (int i = 0; i < Properties.ELITE; i++) {
			logger.trace("Copying individual " + i + " with fitness " + population.get(i).getFitness());
			elite.add((T) population.get(i).clone());
		}
		logger.trace("Done.");
		return elite;
	}

	/**
	 * Create random individuals
	 * 
	 * @return a {@link java.util.List} object.
	 */
	protected List<Chromosome> randomism() {
		logger.debug("Randomism");

		List<Chromosome> randoms = new ArrayList<Chromosome>();

		for (int i = 0; i < Properties.ELITE; i++) {
			randoms.add(chromosomeFactory.getChromosome());
		}
		return randoms;
	}

	/**
	 * update archive fitness functions
	 */
	public void updateFitnessFunctionsAndValues() {
		for (FitnessFunction<T> f : fitnessFunctions) {
			f.updateCoveredGoals();
		}

		// Do we actually have to perform yet another fitness evaluation?
		// Yes, if ARCHIVE has been updated, No otherwise.
		if (!Archive.getArchiveInstance().hasBeenUpdated()) {
			return;
		}

		for (T t : population) {
			for (FitnessFunction<T> fitnessFunction : fitnessFunctions) {
				fitnessFunction.getFitness(t);
			}
		}

		Archive.getArchiveInstance().setHasBeenUpdated(false);
	}

	/**
	 * Penalty if individual is not unique
	 * 
	 * @param individual
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @param generation
	 *            a {@link java.util.List} object.
	 */
	/*
	 * protected void kinCompensation(Chromosome individual, List<Chromosome>
	 * generation) {
	 * 
	 * if (Properties.KINCOMPENSATION >= 1.0) return;
	 * 
	 * boolean unique = true;
	 * 
	 * for (Chromosome other : generation) { if (other == individual) continue;
	 * 
	 * if (other.equals(individual)) { unique = false; break; } }
	 * 
	 * if (!unique) { logger.debug("Applying kin compensation"); if
	 * (fitnessFunction.isMaximizationFunction())
	 * individual.setFitness(individual.getFitness() Properties.KINCOMPENSATION);
	 * else individual.setFitness(individual.getFitness() (2.0 -
	 * Properties.KINCOMPENSATION)); } }
	 */

	/**
	 * Return the individual with the highest fitChromosomeess
	 * 
	 * @return a {@link org.evosuite.ga.Chromosome} object.
	 */
	public T getBestIndividual() {

		if (population.isEmpty()) {
			return this.chromosomeFactory.getChromosome();
		}

		// Assume population is sorted
		return population.get(0);
	}

	/**
	 * Return the individual(s) with the highest fitChromosomeess
	 * 
	 * @return a list of {@link org.evosuite.ga.Chromosome} object(s).
	 */
	public List<T> getBestIndividuals() {

		List<T> bestIndividuals = new ArrayList<T>();

		if (this.population.isEmpty()) {
			bestIndividuals.add(this.chromosomeFactory.getChromosome());
			return bestIndividuals;
		}

		if (Properties.ALGORITHM == Algorithm.NSGAII || Properties.ALGORITHM == Algorithm.SPEA2)
			return population;

		// Assume population is sorted
		bestIndividuals.add(population.get(0));
		return bestIndividuals;
	}

	/**
	 * Write to a file all fitness values of each individuals.
	 *
	 * @param individuals
	 *            a list of {@link org.evosuite.ga.Chromosome} object(s).
	 */
	public void writeIndividuals(List<T> individuals) {
		if (!Properties.WRITE_INDIVIDUALS) {
			return;
		}

		File dir = new File(Properties.REPORT_DIR);
		if (!dir.exists()) {
			if (!dir.mkdirs()) {
				throw new RuntimeException("Cannot create report dir: " + Properties.REPORT_DIR);
			}
		}

		try {
			File populationFile = new File(
					Properties.REPORT_DIR + File.separator + "pareto_" + this.currentIteration + ".csv");
			populationFile.createNewFile();

			FileWriter fw = new FileWriter(populationFile.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			PrintWriter out = new PrintWriter(bw);

			// header
			List<String> l_string = new ArrayList<String>();

			if (Properties.ALGORITHM == Algorithm.NSGAII) {
				l_string.add("rank");
			} else if (Properties.ALGORITHM == Algorithm.SPEA2) {
				l_string.add("strength");
			}

			for (int i = 0; i < this.fitnessFunctions.size(); i++) {
				l_string.add(this.fitnessFunctions.get(i).getClass().getSimpleName());
			}
			out.println(String.join(",", l_string));

			// content
			for (int j = 0; j < individuals.size(); j++) {
				l_string.clear();

				T individual = individuals.get(j);
				if (Properties.ALGORITHM == Algorithm.NSGAII) {
					l_string.add(Integer.toString(individual.getRank()));
				} else if (Properties.ALGORITHM == Algorithm.SPEA2) {
					l_string.add(Double.toString(individual.getDistance()));
				}

				for (int i = 0; i < this.fitnessFunctions.size(); i++) {
					l_string.add(Double.toString(individual.getFitness(this.fitnessFunctions.get(i))));
				}

				out.println(String.join(",", l_string));
			}

			out.close();
			bw.close();
			fw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Set a new factory method
	 * 
	 * @param factory
	 *            a {@link org.evosuite.ga.ChromosomeFactory} object.
	 */
	public void setChromosomeFactory(ChromosomeFactory<T> factory) {
		chromosomeFactory = factory;
	}

	/**
	 * Set a new xover function
	 * 
	 * @param crossover
	 *            a {@link org.evosuite.ga.operators.crossover.CrossOverFunction}
	 *            object.
	 */
	public void setCrossOverFunction(CrossOverFunction crossover) {
		this.crossoverFunction = crossover;
	}

	/**
	 * Add a new search listener
	 * 
	 * @param listener
	 *            a {@link org.evosuite.ga.metaheuristics.SearchListener} object.
	 */
	public void addListener(SearchListener listener) {
		listeners.add(listener);
	}

	/**
	 * Remove a search listener
	 * 
	 * @param listener
	 *            a {@link org.evosuite.ga.metaheuristics.SearchListener} object.
	 */
	public void removeListener(SearchListener listener) {
		listeners.remove(listener);
	}

	/**
	 * Notify all search listeners of search start
	 */
	protected void notifySearchStarted() {
		for (SearchListener listener : listeners) {
			listener.searchStarted(this);
		}
	}

	/**
	 * Notify all search listeners of search end
	 */
	protected void notifySearchFinished() {
		for (SearchListener listener : listeners) {
			listener.searchFinished(this);
		}
	}

	/**
	 * Notify all search listeners of iteration
	 */
	protected void notifyIteration() {
		for (SearchListener listener : listeners) {
			listener.iteration(this);
		}
	}

	/**
	 * Notify all search listeners of fitness evaluation
	 * 
	 * @param chromosome
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 */
	protected void notifyEvaluation(Chromosome chromosome) {
		for (SearchListener listener : listeners) {
			listener.fitnessEvaluation(chromosome);
		}
	}

	/**
	 * Notify all search listeners of a mutation
	 * 
	 * @param chromosome
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 */
	protected void notifyMutation(Chromosome chromosome) {
		for (SearchListener listener : listeners) {
			listener.modification(chromosome);
		}
	}

	/**
	 * Sort the population by fitness
	 * 
	 * WARN: used only with singular objective algorithms, multi-objective
	 * algorithms should implement their own 'sort'
	 */
	protected void sortPopulation() {
		if (Properties.SHUFFLE_GOALS)
			Randomness.shuffle(population);

		if (isMaximizationFunction()) {
			Collections.sort(population, Collections.reverseOrder());
		} else {
			Collections.sort(population);
		}
	}

	/**
	 * Accessor for population Chromosome *
	 * 
	 * @return a {@link java.util.List} object.
	 */
	public List<T> getPopulation() {
		return population;
	}

	/**
	 * Determine if the next generation has reached its size limit
	 * 
	 * @param nextGeneration
	 *            a {@link java.util.List} object.
	 * @return a boolean.
	 */
	public boolean isNextPopulationFull(List<T> nextGeneration) {
		return populationLimit.isPopulationFull(nextGeneration);
	}

	/**
	 * Set a new population limit function
	 * 
	 * @param limit
	 *            a {@link org.evosuite.ga.populationlimit.PopulationLimit} object.
	 */
	public void setPopulationLimit(PopulationLimit limit) {
		this.populationLimit = limit;
	}

	/**
	 * Determine whether any of the stopping conditions hold
	 * 
	 * @return a boolean.
	 */
	public boolean isFinished() {
		for (StoppingCondition c : stoppingConditions) {
			// logger.error(c + " "+ c.getCurrentValue());
			if (c.isFinished())
				return true;
		}
		return false;
	}

	// TODO: Override equals method in StoppingCondition
	/**
	 * <p>
	 * addStoppingCondition
	 * </p>
	 * 
	 * @param condition
	 *            a {@link org.evosuite.ga.stoppingconditions.StoppingCondition}
	 *            object.
	 */
	public void addStoppingCondition(StoppingCondition condition) {
		Iterator<StoppingCondition> it = stoppingConditions.iterator();
		while (it.hasNext()) {
			if (it.next().getClass().equals(condition.getClass())) {
				return;
			}
		}
		logger.debug("Adding new stopping condition");
		stoppingConditions.add(condition);
		addListener(condition);
	}

	public Set<StoppingCondition> getStoppingConditions() {
		return stoppingConditions;
	}

	// TODO: Override equals method in StoppingCondition
	/**
	 * <p>
	 * setStoppingCondition
	 * </p>
	 * 
	 * @param condition
	 *            a {@link org.evosuite.ga.stoppingconditions.StoppingCondition}
	 *            object.
	 */
	public void setStoppingCondition(StoppingCondition condition) {
		stoppingConditions.clear();
		logger.debug("Setting stopping condition");
		stoppingConditions.add(condition);
		addListener(condition);
	}

	/**
	 * <p>
	 * removeStoppingCondition
	 * </p>
	 * 
	 * @param condition
	 *            a {@link org.evosuite.ga.stoppingconditions.StoppingCondition}
	 *            object.
	 */
	public void removeStoppingCondition(StoppingCondition condition) {
		Iterator<StoppingCondition> it = stoppingConditions.iterator();
		while (it.hasNext()) {
			if (it.next().getClass().equals(condition.getClass())) {
				it.remove();
				removeListener(condition);
			}
		}
	}

	/**
	 * <p>
	 * resetStoppingConditions
	 * </p>
	 */
	public void resetStoppingConditions() {
		for (StoppingCondition c : stoppingConditions) {
			c.reset();
		}
	}

	/**
	 * <p>
	 * setStoppingConditionLimit
	 * </p>
	 * 
	 * @param value
	 *            a int.
	 */
	public void setStoppingConditionLimit(int value) {
		for (StoppingCondition c : stoppingConditions) {
			c.setLimit(value);
		}
	}

	protected void updateBestIndividualFromArchive() {
		if (!Properties.TEST_ARCHIVE)
			return;

		T best = Archive.getArchiveInstance().mergeArchiveAndSolution(getBestIndividual());

		// The archive may contain tests evaluated with a fitness function
		// that is not part of the optimization (e.g. ibranch secondary objective)
		Iterator<FitnessFunction<?>> it = best.getCoverageValues().keySet().iterator();
		while (it.hasNext()) {
			FitnessFunction<?> ff = it.next();
			if (!fitnessFunctions.contains(ff))
				it.remove();
		}
		population.add(0, best);
	}

	/**
	 * Returns true if the <code>chromosome1</code> is better or equal than
	 * <code>chromosome2</code> according to the compound fitness function.
	 * 
	 * @param chromosome1
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @param chromosome2
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @return a boolean.
	 */
	protected boolean isBetterOrEqual(Chromosome chromosome1, Chromosome chromosome2) {
		// if (fitnessFunction.isMaximizationFunction()) {
		if (getFitnessFunction().isMaximizationFunction()) {
			return chromosome1.compareTo(chromosome2) >= 0;
		} else {
			return chromosome1.compareTo(chromosome2) <= 0;
		}
	}

	/*
	 * protected boolean isBetter(Chromosome chromosome1, Chromosome chromosome2) {
	 * if (fitnessFunction.isMaximizationFunction()) { return
	 * chromosome1.compareTo(chromosome2) > 0; } else { return
	 * chromosome1.compareTo(chromosome2) < 0; } }
	 */

	/**
	 * <p>
	 * getBest
	 * </p>
	 * 
	 * @param chromosome1
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @param chromosome2
	 *            a {@link org.evosuite.ga.Chromosome} object.
	 * @return a {@link org.evosuite.ga.Chromosome} object.
	 */
	/*
	 * protected Chromosome getBest(Chromosome chromosome1, Chromosome chromosome2)
	 * { if (isBetterOrEqual(chromosome1, chromosome2)) return chromosome1; else
	 * return chromosome2; }
	 */

	/**
	 * Prints out all information regarding this GAs stopping conditions
	 * 
	 * So far only used for testing purposes in TestSuiteGenerator
	 */
	public void printBudget() {
		LoggingUtils.getEvoLogger().info("* GA-Budget:");
		for (StoppingCondition sc : stoppingConditions)
			LoggingUtils.getEvoLogger().info("\t- " + sc.toString());
	}

	/**
	 * <p>
	 * getBudgetString
	 * </p>
	 * 
	 * @return a {@link java.lang.String} object.
	 */
	public String getBudgetString() {
		String r = "";
		for (StoppingCondition sc : stoppingConditions)
			r += sc.toString() + " ";

		return r;
	}

	/**
	 * Returns the progress of the search.
	 * 
	 * @return a value [0.0, 1.0]
	 */
	protected double progress() {
		long totalbudget = 0;
		long currentbudget = 0;

		for (StoppingCondition sc : this.stoppingConditions) {
			if (sc.getLimit() != 0) {
				totalbudget += sc.getLimit();
				currentbudget += sc.getCurrentValue();
			}
		}

		return (double) currentbudget / (double) totalbudget;
	}

	/*
	 * private void writeObject(ObjectOutputStream oos) throws IOException { if
	 * (listeners.contains(SearchStatistics.getInstance())) {
	 * removeListener(SearchStatistics.getInstance()); oos.defaultWriteObject();
	 * oos.writeObject(Boolean.TRUE); // Write/save additional fields
	 * oos.writeObject(SearchStatistics.getInstance()); } else {
	 * oos.defaultWriteObject(); oos.writeObject(Boolean.FALSE); } }
	 * 
	 * // assumes "static java.util.Date aDate;" declared private void
	 * readObject(ObjectInputStream ois) throws ClassNotFoundException, IOException
	 * { ois.defaultReadObject(); listeners = new HashSet<SearchListener>();
	 * stoppingConditions = new HashSet<StoppingCondition>();
	 * 
	 * boolean addStatistics = (Boolean) ois.readObject(); if (addStatistics) {
	 * SearchStatistics.setInstance((SearchStatistics) ois.readObject());
	 * addListener(SearchStatistics.getInstance()); } }
	 */

	//////////////////////////

	// This to generate list of criteria, chosen randomly from the defualt 8
	public Criterion[] getChosenCriteria() {
		Criterion[] originalCriteria = new Criterion[] { Criterion.LINE, Criterion.BRANCH, Criterion.EXCEPTION,
				Criterion.WEAKMUTATION, Criterion.OUTPUT, Criterion.METHOD, Criterion.METHODNOEXCEPTION,
				Criterion.CBRANCH };
		List<Criterion> originalCriteriaList = new ArrayList<Criterion>();
		for (Criterion c : originalCriteria)
			originalCriteriaList.add(c);
		Collections.shuffle(originalCriteriaList);

		/*
		 * random() returns a number between 0.0 and 0.9, 7, becomes 0.0 to 6.93, add 1,
		 * becomes 1.0 to 7.93, the integer will delete the .93.
		 */
		int lngth = (int) (Math.random() * 7 + 1);
		if (lngth == 1)
			System.out.println("h");
		int start = 0;
		if (lngth != 7)
			start = Random.nextInt(7 - lngth);
		System.out.println("lenght" + Integer.toString(lngth));
		System.out.println("lenght" + Integer.toString(start));
		originalCriteriaList = originalCriteriaList.subList(start, lngth);
		Criterion[] subCriteria = new Criterion[lngth];
		originalCriteriaList.toArray(subCriteria);
		return subCriteria;
	}

	public Criterion[] getOneCriteria(int index) {
		// Criterion[] originalCriteria = new Criterion[] {Criterion.LINE,
		// Criterion.BRANCH, Criterion.EXCEPTION, Criterion.WEAKMUTATION,
		// Criterion.OUTPUT, Criterion.METHOD, Criterion.METHODNOEXCEPTION,
		// Criterion.CBRANCH };
		Criterion[] subCriteria = new Criterion[3];
		switch (index) {
		// Combinations of two
		case 0:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE };
			break;
		case 1:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.EXCEPTION };
			break;
		case 2:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.WEAKMUTATION };
			break;
		case 3:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.OUTPUT };
			break;
		case 4:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD };
			break;
		case 5:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.CBRANCH };
			break;
		case 6:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHODNOEXCEPTION };
			break;
		case 7:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.EXCEPTION };
			break;
		case 8:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD };
			break;
		case 9:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.CBRANCH };
			break;
		case 10:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHODNOEXCEPTION };
			break;
		case 11:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.WEAKMUTATION };
			break;
		case 12:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.OUTPUT };
			break;
		case 13:
			subCriteria = new Criterion[] { Criterion.CBRANCH, Criterion.METHOD };
			break;
		case 14:
			subCriteria = new Criterion[] { Criterion.CBRANCH, Criterion.EXCEPTION };
			break;
		case 15:
			subCriteria = new Criterion[] { Criterion.CBRANCH, Criterion.METHODNOEXCEPTION };
			break;
		case 16:
			subCriteria = new Criterion[] { Criterion.CBRANCH, Criterion.WEAKMUTATION };
			break;
		case 17:
			subCriteria = new Criterion[] { Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 18:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.EXCEPTION };
			break;
		case 19:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.METHODNOEXCEPTION };
			break;
		case 20:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.WEAKMUTATION };
			break;
		case 21:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 22:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 23:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION };
			break;
		case 24:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 25:
			subCriteria = new Criterion[] { Criterion.METHODNOEXCEPTION, Criterion.WEAKMUTATION };
			break;
		case 26:
			subCriteria = new Criterion[] { Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 27:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.OUTPUT };
			break;

		// Combinations of three

		case 28:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.OUTPUT };
			break;
		case 29:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.CBRANCH };
			break;
		case 30:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.METHODNOEXCEPTION };
			break;
		case 31:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.EXCEPTION };
			break;
		case 32:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.WEAKMUTATION };
			break;
		case 33:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.LINE, Criterion.METHOD };
			break;
		case 34:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 35:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD, Criterion.CBRANCH };
			break;
		case 36:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD, Criterion.METHODNOEXCEPTION };
			break;
		case 37:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD, Criterion.EXCEPTION };
			break;
		case 38:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHOD, Criterion.WEAKMUTATION };
			break;
		case 39:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.WEAKMUTATION, Criterion.OUTPUT };
			break;
		case 40:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.WEAKMUTATION, Criterion.CBRANCH };
			break;
		case 41:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION };
			break;
		case 42:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.WEAKMUTATION, Criterion.EXCEPTION };
			break;
		case 43:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 44:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.EXCEPTION, Criterion.CBRANCH };
			break;
		case 45:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 46:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 47:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH };
			break;
		case 48:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 49:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 50:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD, Criterion.CBRANCH };
			break;
		case 51:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD, Criterion.METHODNOEXCEPTION };
			break;
		case 52:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD, Criterion.EXCEPTION };
			break;
		case 53:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHOD, Criterion.WEAKMUTATION };
			break;
		case 54:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.WEAKMUTATION, Criterion.OUTPUT };
			break;
		case 55:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.WEAKMUTATION, Criterion.CBRANCH };
			break;
		case 56:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION };
			break;
		case 57:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.WEAKMUTATION, Criterion.EXCEPTION };
			break;
		case 58:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 59:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.EXCEPTION, Criterion.CBRANCH };
			break;
		case 60:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 61:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 62:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH };
			break;
		case 63:
			subCriteria = new Criterion[] { Criterion.LINE, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 64:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.WEAKMUTATION, Criterion.OUTPUT };
			break;
		case 65:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.WEAKMUTATION, Criterion.CBRANCH };
			break;
		case 66:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION };
			break;
		case 67:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.WEAKMUTATION, Criterion.EXCEPTION };
			break;
		case 68:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 69:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.EXCEPTION, Criterion.CBRANCH };
			break;
		case 70:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 71:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 72:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH };
			break;
		case 73:
			subCriteria = new Criterion[] { Criterion.METHOD, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 74:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 75:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.EXCEPTION, Criterion.CBRANCH };
			break;
		case 76:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 77:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 78:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH };
			break;
		case 79:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 80:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 81:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH };
			break;
		case 82:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 83:
			subCriteria = new Criterion[] { Criterion.METHODNOEXCEPTION, Criterion.CBRANCH, Criterion.OUTPUT };
			break;

		// one
		case 84:
			subCriteria = new Criterion[] { Criterion.BRANCH };
			break;
		case 85:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION };
			break;
		case 86:
			subCriteria = new Criterion[] { Criterion.EXCEPTION };
			break;
		case 87:
			subCriteria = new Criterion[] { Criterion.METHODNOEXCEPTION };
			break;
		case 88:
			subCriteria = new Criterion[] { Criterion.OUTPUT };
			break;
		case 89:
			subCriteria = new Criterion[] { Criterion.CBRANCH };
			break;
		case 90:
			subCriteria = new Criterion[] { Criterion.METHOD };
			break;
		case 91:
			subCriteria = new Criterion[] { Criterion.LINE };
			break;

		}

		// Criterion[] subCriteria = new Criterion[] {originalCriteria[index]};
		return subCriteria;
	}

	public Criterion[] getOneCriteriaWithException(int index) {
		// Criterion[] originalCriteria = new Criterion[] {Criterion.LINE,
		// Criterion.BRANCH, Criterion.EXCEPTION, Criterion.WEAKMUTATION,
		// Criterion.OUTPUT, Criterion.METHOD, Criterion.METHODNOEXCEPTION,
		// Criterion.CBRANCH };
		Criterion[] subCriteria = new Criterion[3];
		switch (index) {
		case 0:
			subCriteria = new Criterion[] { Criterion.BRANCH, Criterion.EXCEPTION, Criterion.LINE, Criterion.CBRANCH };
			break;
		case 1:
			subCriteria = new Criterion[] { Criterion.WEAKMUTATION, Criterion.EXCEPTION, Criterion.LINE };
			break;
		case 2:
			subCriteria = new Criterion[] { Criterion.EXCEPTION };
			break;
		case 3:
			subCriteria = new Criterion[] { Criterion.METHODNOEXCEPTION };
			break;
		case 4:
			subCriteria = new Criterion[] { Criterion.OUTPUT, Criterion.EXCEPTION, Criterion.LINE };
			break;
		case 5:
			subCriteria = new Criterion[] { Criterion.CBRANCH };
			break;
		case 6:
			subCriteria = new Criterion[] { Criterion.METHOD };
			break;
		case 7:
			subCriteria = new Criterion[] { Criterion.LINE };
			break;

		case 8:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH };
			break;
		case 9:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION };
			break;
		case 10:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION };
			break;
		case 11:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.OUTPUT };
			break;
		case 12:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH };
			break;
		case 13:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD };
			break;
		case 14:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE };
			break;

		case 15:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE };
			break;
		case 16:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.WEAKMUTATION };
			break;
		case 17:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.OUTPUT };
			break;
		case 18:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD };
			break;
		case 19:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.CBRANCH };
			break;
		case 20:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHODNOEXCEPTION };
			break;
		case 21:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHOD };
			break;
		case 22:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.CBRANCH };
			break;
		case 23:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHODNOEXCEPTION };
			break;
		case 24:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.WEAKMUTATION };
			break;
		case 25:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.OUTPUT };
			break;
		case 26:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH, Criterion.METHOD };
			break;
		case 27:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH, Criterion.METHODNOEXCEPTION };
			break;
		case 28:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH, Criterion.WEAKMUTATION };
			break;
		case 29:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 30:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.METHODNOEXCEPTION };
			break;
		case 31:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.WEAKMUTATION };
			break;
		case 32:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 33:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION, Criterion.WEAKMUTATION };
			break;
		case 34:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION, Criterion.OUTPUT };
			break;
		case 35:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION, Criterion.OUTPUT };
			break;

		// Combinations of three
		case 36:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE, Criterion.OUTPUT };
			break;
		case 37:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE, Criterion.CBRANCH };
			break;
		case 38:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE,
					Criterion.METHODNOEXCEPTION };
			break;
		case 39:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE,
					Criterion.WEAKMUTATION };
			break;
		case 40:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.LINE, Criterion.METHOD };
			break;
		case 41:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 42:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD,
					Criterion.CBRANCH };
			break;
		case 43:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD,
					Criterion.METHODNOEXCEPTION };
			break;
		case 44:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD,
					Criterion.EXCEPTION };
			break;
		case 45:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHOD,
					Criterion.WEAKMUTATION };
			break;
		case 46:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.WEAKMUTATION,
					Criterion.OUTPUT };
			break;
		case 47:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.WEAKMUTATION,
					Criterion.CBRANCH };
			break;
		case 48:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.WEAKMUTATION,
					Criterion.METHODNOEXCEPTION };
			break;
		case 49:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHODNOEXCEPTION,
					Criterion.OUTPUT };
			break;
		case 50:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.METHODNOEXCEPTION,
					Criterion.CBRANCH };
			break;
		case 51:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.BRANCH, Criterion.CBRANCH,
					Criterion.OUTPUT };
			break;
		case 52:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHOD, Criterion.OUTPUT };
			break;
		case 53:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHOD, Criterion.CBRANCH };
			break;
		case 54:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHOD,
					Criterion.METHODNOEXCEPTION };
			break;
		case 55:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHOD,
					Criterion.WEAKMUTATION };
			break;
		case 56:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.WEAKMUTATION,
					Criterion.OUTPUT };
			break;
		case 57:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.WEAKMUTATION,
					Criterion.CBRANCH };
			break;
		case 58:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.WEAKMUTATION,
					Criterion.METHODNOEXCEPTION };
			break;
		case 59:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHODNOEXCEPTION,
					Criterion.OUTPUT };
			break;
		case 60:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.METHODNOEXCEPTION,
					Criterion.CBRANCH };
			break;
		case 61:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.LINE, Criterion.CBRANCH, Criterion.OUTPUT };
			break;
		case 62:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.WEAKMUTATION,
					Criterion.OUTPUT };
			break;
		case 63:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.WEAKMUTATION,
					Criterion.CBRANCH };
			break;
		case 64:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.WEAKMUTATION,
					Criterion.METHODNOEXCEPTION };
			break;
		case 65:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.METHODNOEXCEPTION,
					Criterion.OUTPUT };
			break;
		case 66:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.METHODNOEXCEPTION,
					Criterion.CBRANCH };
			break;
		case 67:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHOD, Criterion.CBRANCH,
					Criterion.OUTPUT };
			break;
		case 68:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION,
					Criterion.OUTPUT };
			break;
		case 69:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION, Criterion.METHODNOEXCEPTION,
					Criterion.CBRANCH };
			break;
		case 70:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.WEAKMUTATION, Criterion.CBRANCH,
					Criterion.OUTPUT };
			break;
		case 71:
			subCriteria = new Criterion[] { Criterion.EXCEPTION, Criterion.METHODNOEXCEPTION, Criterion.CBRANCH,
					Criterion.OUTPUT };
			break;

		}

		// Criterion[] subCriteria = new Criterion[] {originalCriteria[index]};
		return subCriteria;
	}

	public void removeFitnessFunction(Criterion[] temp) {
		// TestsArchive.instance.reset();
		Properties.CRITERION = temp;
		List<TestSuiteFitnessFunction> ffs = getFitnessFunctions1();
		this.fitnessFunctions.clear();
		this.addFitnessFunctions((List) ffs);
		List<TestFitnessFunction> goals = getGoals(true);
		// deleteFitnessFromPopulation();
		calculateFitnessAndSortPopulation();
		// updateFitnessFunctionsAndValues(true);
		deleteUpdatFitnessFromPopulation();
	}

	/// this to remove the deleted fitness functions from the population
	protected void deleteFitnessFromPopulation() {
		for (T individual : population) {
			individual.clearFitnessVlaue();
			for (FitnessFunction<?> fitnessFunction : this.fitnessFunctions) {
				individual.addFitness(fitnessFunction);
			}
		}
	}

	/**
	 * update archive fitness functions
	 */
	protected void updateFitnessFunctionsAndValues(boolean fitnessNeedsUpdating) {
		for (FitnessFunction<T> f : fitnessFunctions) {
			f.updateCoveredGoals();
		}
		for (T t : population) {
			for (FitnessFunction<T> fitnessFunction : fitnessFunctions) {
				fitnessFunction.getFitness(t);
			}
		}
	}

	protected void deleteUpdatFitnessFromPopulation() {
		for (T individual : population) {
			individual.clearFitnessVlaue();
			for (FitnessFunction<T> fitnessFunction : this.fitnessFunctions) {
				individual.addFitness(fitnessFunction);
				fitnessFunction.updateCoveredGoals();
				fitnessFunction.getFitness(individual);
			}
		}
	}

	/**
	 * 
	 * 
	 * /** Convert criterion names to test suite fitness functions
	 * 
	 * @return
	 */
	protected List<TestSuiteFitnessFunction> getFitnessFunctions1() {
		List<TestSuiteFitnessFunction> ffs = new ArrayList<TestSuiteFitnessFunction>();
		for (int i = 0; i < Properties.CRITERION.length; i++) {
			TestSuiteFitnessFunction newFunction = FitnessFunctions.getFitnessFunction(Properties.CRITERION[i]);

			// If this is compositional fitness, we need to make sure
			// that all functions are consistently minimization or
			// maximization functions
			if (Properties.ALGORITHM != Algorithm.NSGAII && Properties.ALGORITHM != Algorithm.SPEA2) {
				for (TestSuiteFitnessFunction oldFunction : ffs) {
					if (oldFunction.isMaximizationFunction() != newFunction.isMaximizationFunction()) {
						StringBuffer sb = new StringBuffer();
						sb.append("* Invalid combination of fitness functions: ");
						sb.append(oldFunction.toString());
						if (oldFunction.isMaximizationFunction())
							sb.append(" is a maximization function ");
						else
							sb.append(" is a minimization function ");
						sb.append(" but ");
						sb.append(newFunction.toString());
						if (newFunction.isMaximizationFunction())
							sb.append(" is a maximization function ");
						else
							sb.append(" is a minimization function ");
						LoggingUtils.getEvoLogger().info(sb.toString());
						throw new RuntimeException("Invalid combination of fitness functions");
					}
				}
			}
			ffs.add(newFunction);

		}

		return ffs;
	}

	/**
	 * Convert criterion names to factories for test case fitness functions
	 * 
	 * @return
	 */
	public static List<TestFitnessFactory<? extends TestFitnessFunction>> getFitnessFactories() {
		List<TestFitnessFactory<? extends TestFitnessFunction>> goalsFactory = new ArrayList<TestFitnessFactory<? extends TestFitnessFunction>>();
		for (int i = 0; i < Properties.CRITERION.length; i++) {
			goalsFactory.add(FitnessFunctions.getFitnessFactory(Properties.CRITERION[i]));
		}

		return goalsFactory;
	}

	private List<TestFitnessFunction> getGoals(boolean verbose) {
		List<TestFitnessFactory<? extends TestFitnessFunction>> goalFactories = getFitnessFactories();
		List<TestFitnessFunction> goals = new ArrayList<>();

		if (goalFactories.size() == 1) {
			TestFitnessFactory<? extends TestFitnessFunction> factory = goalFactories.iterator().next();
			goals.addAll(factory.getCoverageGoals());

			if (verbose) {
				LoggingUtils.getEvoLogger().info("* Total number of test goals: {}", factory.getCoverageGoals().size());
				if (Properties.PRINT_GOALS) {
					for (TestFitnessFunction goal : factory.getCoverageGoals())
						LoggingUtils.getEvoLogger().info("" + goal.toString());
				}
			}
		} else {
			if (verbose) {
				LoggingUtils.getEvoLogger().info("* Total number of test goals: ");
			}

			for (TestFitnessFactory<? extends TestFitnessFunction> goalFactory : goalFactories) {
				goals.addAll(goalFactory.getCoverageGoals());

				if (verbose) {
					LoggingUtils.getEvoLogger()
							.info("  - " + goalFactory.getClass().getSimpleName().replace("CoverageFactory", "") + " "
									+ goalFactory.getCoverageGoals().size());
					if (Properties.PRINT_GOALS) {
						for (TestFitnessFunction goal : goalFactory.getCoverageGoals())
							LoggingUtils.getEvoLogger().info("" + goal.toString());
					}
				}
			}
		}
		return null;
	}

	protected void sortPopulationNew() {
		if (Properties.SHUFFLE_GOALS)
			Randomness.shuffle(population);

		if (isMaximizationFunction()) {
			Collections.sort(population, Collections.reverseOrder());
		} else {
			T tMax = population.get(0);
			TestSuiteChromosome max = (TestSuiteChromosome) population.get(0);
			int[] fsmax = CoverageCriteriaAnalyzer.analyzeCoverageNew(max, Criterion.EXCEPTION);
			int fmax = fsmax[0] + fsmax[1];

			for (int i = 1; i < population.size(); i++) {

				TestSuiteChromosome tmp = (TestSuiteChromosome) population.get(i);
				int[] fs = CoverageCriteriaAnalyzer.analyzeCoverageNew(tmp, Criterion.EXCEPTION);
				int f = fs[0] + fs[1];

				if (f > fmax) {
					tMax = population.get(i);
					max = tmp;
					fsmax = fs;
					fmax = f;
				}
			}
			T tmp1 = population.get(0);
			int ind = population.indexOf(tMax);
			population.set(0, tMax);
			population.set(ind, tmp1);
		}
	}

	///////////////////////

}
