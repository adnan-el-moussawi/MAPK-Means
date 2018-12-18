package clustering;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import comparator.MinkowskiTab;
import utilities.ArraysUtilities;



/**
 * MAPK-Means : Metric Learning Preferential K-Means (El Moussawi et al., 2016)
 * <br>This class implements MAPK-Means clustering algorithm, an algorithm 
 * with quantitative user preferences on attributes (See the abstract in the following).
 * <br>
 * <br><b>Abstract:</b>
 * <br>The paper proposes a new semi-supervised clustering framework to represent 
 * and integrate quantitative preferences on attributes. A new metric learning 
 * algorithm is derived that achieves a compromise clustering between a data-driven 
 * and a user-driven solution and converges with a good complexity. We observe 
 * experimentally that the addition of preferences may be essential to achieve 
 * a better clustering. We also show that our approach performs better than 
 * the state-of-the art algorithms.<br>
 * <br>We encourage you to cite our paper
 * <br>Cite the paper as: (APA): El Moussawi, A., Cheriat, A., Giacometti, A., Labroche, N., 
 * & Soulet, A. (2016, November). Clustering with Quantitative User Preferences 
 * on Attributes. In Tools with Artificial Intelligence (ICTAI), 2016 IEEE 28th 
 * International Conference on (pp. 383-387). IEEE.
 * <br>Cite (BibTex):
 * <p>
 * \@inproceedings{el2016clustering,
 * 	title={Clustering with Quantitative User Preferences on Attributes},
 * 	author={El Moussawi, Adnan and Cheriat, Ahmed and Giacometti, Arnaud and Labroche, Nicolas and Soulet, Arnaud},
 * 	booktitle={Tools with Artificial Intelligence (ICTAI), 2016 IEEE 28th International Conference on},
 * 	pages={383--387},
 * 	year={2016},
 * 	organization={IEEE}
 * }
 * </p>
 * <br>http://doi.ieeecomputersociety.org/10.1109/ICTAI.2016.0065
 * <br>In the following, the objective function of the algorithm :
 * <br>I = alpha*Z*Intra + (1 - alpha)[ kappa*KL(W*,W) + (1 - kappa)*KL(U|W) ) + lambda*(sum(w_i) - 1]
 * <p>
 * This function is a compromise between : (i) the weighted intra-cluster distance, 
 * (ii) the KL divergence between the learned metric and the user preferences 
 * (iii) and the KL divergence between the user preferences and the learned metric 
 * </p>
 * <p>
 * The parameter alpha is the weight (in [0, 1]) of the intra-distance, with 0.5 the default value
 * (proved experimentally, see the paper). The parameter Z is a normalization factor to balance 
 * the inta-cluster with other terms. Finally, kappa is a weight to express the user confidence level 
 * in his or her preferences, its can ranges between 0 and 1. The higher kappa is, more the preferences 
 * are considered to build the clustering. 
 * </p>
 * @author Adnan EL MOUSSAWI
 *
 */
public class MAPKMeans {
	/**
	 * precision
	 */
	public static final double epsilon = 0.000001;
	
	/** 
	 * MIN difference of objective function values for convergence
	 */
	protected final double m_ObjFunConvergenceDifference = 0.000001;
	
	/**
	 * used to save the value of the objective function at each iteration 
	 */
	protected double m_CurrentObjective = Double.MAX_VALUE;
	
	/**
	 * used to save the value of the objective function at the previous iteration
	 */
	protected double m_OldObjective;
	
	
	/**
	 * number of desired cluster
	 */
	private int k;
	/**
	 * the data partition
	 */
	private int[] partition;
	
	/**
	 * vector of clusters size
	 */
	private int[] clustersEffective;
	
	/**
	 * maximal allowed number of iteration
	 */
	private int nbitermax=1000;
	/**
	 * number of iteration to the convergence
	 */
	private int nbiter;
	
	/**
	 * learning weight
	 */
	private double alpha=0.5;
	/**
	 * the normalization factor between the distance and the learning part
	 */
	private double normalizationFactor;
	/**
	 * confidence value
	 */
	private double confidenceLevel=1;
	
	/**
	 * Lagrangian value
	 */
	private double lagrangian;
	
	/**
	 * the user preferences vector
	 */
	private double[] preferencesVector;
	
	/**
	 * the learning vector
	 */
	private double[] learningVector;
	
	/**
	 * the vector of the uniform distribution of weights
	 */
	private double[] uniformDistribution;
	
	/**
	 * number of the dimensions
	 */
	private int nbDimensions;
	
	/**
	 * the data set 
	 */
	private List<double[]> data;
	
	/**
	 * number of data points
	 */
	private int nbData;
	

	/**
	 * the total weighted intra-cluster distance
	 */
	private double intraClusterDistance;
	
	/** 
	 * the Kulback-Leibler divergence between the preferences and the learned weights
	 */
	private double divergenceWithPreferences;
	/**
	 * the Kulback-Leibler divergence between the preferences and the uniform distribution
	 */
	private double divergenceWithUniformDistribution;
	
	
	
	/*private List<Double> objectiveFunctionBeforeAssignment = new ArrayList<Double>(),
	objectiveFBeforeUpdate = new ArrayList<Double>(),
	objectiveFAfterNormalisation = new ArrayList<Double>();*/
	
	/**
	 * the computed centers
	 */
	private List<double[]> centers;

	/**
	 * the initialization class {@link KMeansCentroidInitialization}
	 */
	private KMeansCentroidInitialization initialization = new KMeansCentroidInitialization();
	/**
	 * the data comparator {@link Comparator} to use the distance functions
	 */
	private MinkowskiTab cp = new MinkowskiTab(2);
	
	/**
	 * Main clustering function, the learned vector is initialized with the user preferences, 
	 * and initial centers are computed with Weighted KMeans++
	 * @param data : data set
	 * @param kint number of desired cluster
	 * @param iteration maximal number of iteration until convergence
	 * @param weights preferences vector of weights
	 * @param alpha learning term weight
	 * @param kappa user confidence weight
	 */
	public List<Double> cluster(List<double[]> data, int kint, int iteration,
			double[] weights,	double alpha, double kappa){
		
		double[] tmp = new double[weights.length];
		for(int i =0; i<weights.length; i++){
			tmp[i] = kappa*weights[i] + (1.0-kappa)/weights.length; 
		}
		List<double[]> centers = initialization.initializationKMeansPlusPlus(data, kint, tmp);
		//List<double[]> centers = initialization.AnovaPartitioningInitialization(data, kint);
		
		initializeVariables(data, centers, kint, iteration, weights, alpha, kappa);		
		return clusteringStart();
	}
	/*
	public List<Double> cluster2(List<double[]> data, int kint, int iteration,
			double[] weights,	double alpha, double kappa){
		double[] tmp = new double[weights.length];
		for(int i =0; i<weights.length; i++){
			tmp[i] = kappa*weights[i] + (1.0-kappa)/weights.length; 
		}
		List<double[]> centers = initialization.initializationKMeansPlusPlus(data, kint, tmp);
		//List<double[]> centers = initialization.AnovaPartitioningInitialization(data, kint);
		kint = centers.size();
		return cluster(data, centers, kint, iteration, weights, alpha, kappa);
	}*/
	
	/**
	 * Main clustering function, the learned vector is initialized with the user preferences
	 * @param data : data set
	 * @param centers initial centers
	 * @param kint number of desired cluster
	 * @param iteration maximal number of iteration until convergence, the default value 100
	 * @param weights preferences vector of weights
	 * @param alpha learning term weight, the default value is 0.5
	 * @param kappa user confidence weight
	 */
	public List<Double> cluster(List<double[]> data, List<double[]> centers, int kint, int iteration,
			double[] weights,double alpha, double kappa){
		
		initializeVariables(data, centers, kint, iteration, weights, alpha, kappa);
		
		return clusteringStart();
	}

	/**
	 * Main clustering function, the learned vector is initialized with the uniform distribution
	 * @param data : data set
	 * @param centers initial centers
	 * @param kint number of desired cluster
	 * @param iteration maximal number of iteration until convergence
	 * @param weights preferences vector of weights
	 * @param alpha learning term weight
	 * @param kappa user confidence weight
	 */
	public List<Double> clusterUniformInit(List<double[]> data, List<double[]> centers, int kint, int iteration,
			double[] weights,	double alpha, double kappa){
		/**
		 * algorithm initialization
		 */
		initializeVariables(data, centers, kint, iteration, weights, alpha, kappa);
		learningVector = uniformDistribution.clone();
		return clusteringStart();
	}

	/**
	 * Main clustering function, the learned vector is initialized with the uniform distribution
	 * @param data : data set
	 * @param kint number of desired cluster
	 * @param iteration maximal number of iteration until convergence
	 * @param weights preferences vector of weights
	 * @param alpha learning term weight
	 * @param kappa user confidence weight
	 */
	public List<Double> clusterUniformInit(List<double[]> data, int kint, int iteration,
			double[] weights,	double alpha, double kappa){
		/**
		 * algorithm initialization
		 */
		List<double[]> centers = initialization.initializationKMeansPlusPlus(data, kint);
		
		initializeVariables(data, centers, kint, iteration, weights, alpha, kappa);
		learningVector = uniformDistribution.clone();
		return clusteringStart();
	}
	
	/**
	 * Initialization
	 * @param data : data set
	 * @param centers initial centers
	 * @param cp : comparator {@link Comparator}
	 * @param kint number of desired cluster
	 * @param iteration maximal number of iteration until convergence
	 * @param weights preferences vector of weights
	 * @param alpha learning term weight
	 * @param kappa user confidence weight
	 */
	private void initializeVariables(List<double[]> data, List<double[]> centers, int kint, int iteration,
			double[] weights,	double alpha, double kappa) {
		// TODO Auto-generated method stub
		this.k = kint;
		this.nbitermax = iteration;
		
		this.data = data;
		this.nbData = data.size();
		this.nbDimensions = data.get(0).length;
		
		this.alpha = alpha;
		this.confidenceLevel = kappa;
		this.centers = centers;
		
		this.uniformDistribution = new double[nbDimensions];
		
		Arrays.fill(uniformDistribution, 1.0/nbDimensions);
		 
		this.learningVector =  new double[nbDimensions];
		if(weights == null) {
			this.preferencesVector = uniformDistribution.clone();
			this.learningVector = uniformDistribution.clone();
		} else {
			this.preferencesVector = weights;
			for(int i =0; i<nbDimensions; i++){
				learningVector[i] = kappa*preferencesVector[i] + (1.0-kappa)/nbDimensions; 
			}
			learningVector = ArraysUtilities.normalizeArray(learningVector);
		}
		
		
		partition = new int[nbData];
		nbiter = 0;
	}

	/**
	 * Main clustering function, iterate until convergence
	 * @return The values of the objective function
	 */
	private List<Double> clusteringStart(){
		
		boolean converged = false;
		double dist;
		int best;
		
		List<Double> objectiveFunction = new ArrayList<Double>();
		
		double numerator_sup = ((double)nbDimensions)*(1.0-alpha)*
				( confidenceLevel*ArraysUtilities.getMaxArray(preferencesVector) + (1.0 - confidenceLevel) );
		while (!converged && nbiter < nbitermax){
			clustersEffective = new int[k];
			double tmp;
			/**
			 * assignment step
			 */
			for (int i = 0; i < nbData; i++){
				dist = cp.weightedDistance2(data.get(i), centers.get(0), learningVector);
				
				best = 0;
				for (int c = 1; c < k; c++){
					tmp = cp.weightedDistance2(data.get(i), centers.get(c), learningVector);
					if (tmp < dist) {
						dist = tmp; 
						best = c;
					}
				}
				partition[i] = best;
				clustersEffective[best]++;
				//break;
			}
			
			/**
			 * centers update
			 */
			centers = new ArrayList<double[]>();
			for (int i = 0; i < k; i++)
				centers.add(new double[data.get(0).length]);
			for (int i = 0; i < nbData; i++)
				add(centers, data.get(i), partition[i]);
			for (int i = 0; i < k; i++)
				div(centers, clustersEffective[i], i);
			
			double[] s = computeIntraClusterVariance();
			
			/*
			 * set the normalization factor
			 */
			if (nbiter==0) {
				setNormalizationFactor(s);
			}
			
			/*
			 * metric learning step
			 */
			if(this.alpha!=1){

				double denominator_inf = alpha*normalizationFactor*ArraysUtilities.getMinArray(s);
				double lagrang_inf = (- denominator_inf);
				double lagrang_sup =  numerator_sup  - denominator_inf;
				lagrangian = (double) computeLagrangeValue(s, lagrang_inf, lagrang_sup);				
				learningVector = updateLearningVector(s);
			}
			
			/*
			 * update the objective function and check the convergence
			 */
			m_OldObjective = m_CurrentObjective;
			m_CurrentObjective = updateObjectiveFunction();
			converged = convergenceCheck(m_CurrentObjective, m_OldObjective);
			
			//System.out.println("objective function : "+m_CurrentObjective);
			
			objectiveFunction.add(new Double(getObjectiveFunction()));
			/* Increment the number of current iteration */
			nbiter++;			
		}
		intraClusterDistance = computeTotalWeightedIntraClusterDistance();
		divergenceWithPreferences = KullbackLeibler(preferencesVector, learningVector);
		divergenceWithUniformDistribution = KullbackLeibler(uniformDistribution, learningVector);
		return objectiveFunction;
	}
	
	/**
	 * compute the normalization factor
	 * @param s
	 */
	private void setNormalizationFactor(double[] s) {
		// TODO Auto-generated method stub
		double tmp = 0.0;
		for(int i=0; i<s.length; i++)
			tmp += (confidenceLevel*preferencesVector[i]+(1-confidenceLevel)/(double) nbDimensions)/s[i];
		this.normalizationFactor = tmp;
	}

	/**
	 * Compute the euclidean intra-cluster variance for the clustering on each dimension
	 * @return vector of clustering variance on each dimension
	 */
	private double[] computeIntraClusterVariance() {
		// TODO Auto-generated method stub
		
		double[] s = new double[nbDimensions];
		for(int d=0; d<nbDimensions; d++){
			s[d]=0.0;
			for(int i=0; i<nbData; i++){
				for(int c=0; c<centers.size(); c++){
					if(partition[i] == c){
						s[d] += Math.pow((data.get(i)[d] - centers.get(c)[d]), 2.0);
					}
				}
			}
			//System.out.println("S["+d+"] =  "+s[d]);
		}
		return s.clone();
	}
	
	
	/**
	 * Update the objective function
	 * @return the new value of the objective function
	 */
	private double updateObjectiveFunction(){		
		double tmp= 0.0;
		
		for(int i=0; i<nbData; i++){
			for(int c=0; c<centers.size(); c++){
				if(partition[i]==c){
					tmp += cp.weightedDistance2(data.get(i), centers.get(c), learningVector);
				}
			}
		}
		
		tmp = alpha*normalizationFactor*tmp + (1.0 - alpha)*confidenceLevel*KullbackLeibler(preferencesVector, learningVector)
				+ (1.0 - alpha)*(1.0 - confidenceLevel)*KullbackLeibler(uniformDistribution, learningVector)
				+ lagrangian*(ArraysUtilities.getSumArray(learningVector)-1.0);
		
		return tmp;
	}
	
	/**
	 * Compute the Lagrangian using  a dichotomic search
	 * @param s vector of the intra-cluster euclidean variance on each dimension
	 * @param lagrang_inf lagrangian Inf born
	 * @param lagrang_sup lagrangian Sup born
	 * @return the Lagrangian value
	 */
    
	 private double computeLagrangeValue(double[] s, double lagrang_inf, double lagrang_sup) {
		double lambda_0_plus = 1E-17;
		
		double temp_lambda = (lagrang_inf + lambda_0_plus + lagrang_sup)/2.0;
		double sum_ai = 0.0;		
		
		//false while value not founded
		boolean founded = false, precisionAttend = false;
		
		while(!founded && !precisionAttend){
			sum_ai = 0.0;
			for(int d=0; d<nbDimensions; d++){
				sum_ai += (1.0-alpha)*( confidenceLevel*preferencesVector[d] + (1.0 - confidenceLevel)/( (double) nbDimensions ) )
						/( alpha*normalizationFactor*s[d] + temp_lambda );	
			}

			founded = (Math.abs(sum_ai - 1.0) <= 1E-15);
			if(founded) return temp_lambda;
			
			//System.out.println("sum_ai : "+sum_ai);
			//System.out.println(temp_lambda);
			
			if (sum_ai < 1.0){
				lagrang_sup = temp_lambda;
			}
			else if( sum_ai > 1.0){
				lagrang_inf = temp_lambda;
			}
			else{
				return temp_lambda;
			}
			
			precisionAttend = Math.abs(lagrang_sup - lagrang_inf) <= (1.0E-15);
			
			temp_lambda = (lagrang_inf + lagrang_sup)/2.0;
		}
		
		//System.out.println(temp_lambda);
		
		return (double) temp_lambda;
	}
	 
	
	/**
	 * Update the metric (the vector of learned weights)
	 * @param s The vector of the intra-cluster euclidean variance on each dimension
	 * @return
	 */
	private double[] updateLearningVector(double[] s) {
		// TODO Auto-generated method stub
		double[] a = new double[nbDimensions];
		for(int i=0; i<nbDimensions; i++){
			if(alpha != 1.0){
				a[i] = (1.0-alpha)*( confidenceLevel*preferencesVector[i] + (1.0 - confidenceLevel)/( (double) nbDimensions ) )
						/( alpha*normalizationFactor*s[i] + lagrangian );
			}
			else{
				a[i] = 1.0/s[i];
			}
		}
	
		//System.out.println(Arrays.toString(a));
		return a.clone();
	}

	/**
	 * Compute the Kullback-Leibler divergence between two vectors
	 * @param a1 
	 * @param a2 
	 * @return The KL divergence between two vectors
	 */
	private double KullbackLeibler (double[] a1, double[] a2){
		//System.out.println(Arrays.toString(a1));
		//System.out.println(Arrays.toString(a2));
		double tmp = 0;
		for(int i=0; i<a1.length; i++){
			if(a1[i]==0.0) continue;//a*.log(a*/a_i)=0
			else{
				tmp += a1[i]*Math.log(a1[i]/a2[i]);
			}
		}
		return tmp;
	}
	
	/**
	 * Check the convergence between two values
	 * @param m_CurrentObjective 
	 * @param m_OldObjective 
	 * @return true if the values are closed (with precision equal to 10^{-6})
	 */
	private boolean convergenceCheck(double m_CurrentObjective, double m_OldObjective) {
		// TODO Auto-generated method stub
		boolean converged = false;
	
		// Convergence check
		if(Math.abs(m_OldObjective - m_CurrentObjective) < m_ObjFunConvergenceDifference) {
			//System.out.println("Final objective function is: " + m_CurrentObjective);
			converged = true;
		}
	
		return converged;
	}
	
	/**
	 * Compute the total weighted intra-cluster variance
	 * @return
	 */
	private double computeTotalWeightedIntraClusterDistance(){
		double tmp=0.0;
		for(int i=0; i<nbData; i++){
			for(int c=0; c<centers.size(); c++){
				if(partition[i]==c){
					tmp += cp.weightedDistance2(data.get(i), centers.get(c), learningVector);
				}
			}			
		}
		return tmp;
	}
	
	/**
	 * Get the total weighted intra-cluster variance
	 * @return
	 */
	public double getIntraClusterDistance(){
		return intraClusterDistance;
	}

	/**
	 * Get the divergence between the learned vector and the user preferences
	 * <br>KL(W||W^*)
	 * @return
	 */
	public double getDivergenceWithPreferences(){
		return this.divergenceWithPreferences;
	}
	
	/**
	 * Get the divergence between the learned vector and the uniform distribution
	 * <br>KL(W||U)
	 * @return
	 */
	public double getDivergenceWithUniformDistribution(){
		return this.divergenceWithUniformDistribution;
	}
	/**
	 * Get the first part of objective function : 
	 * <br>(\alpha \times \zFactor \times weightedIntra)
	 * @return
	 */
	public double getDistanceTerm(){
		return alpha*normalizationFactor*intraClusterDistance;
	}
	
	/**
	 * Get the second part (learning part) of the objective function
	 * <br>(\kappa \times KL(W||W^*) + (1-\kappa) \times KL(W||U))
	 * @return
	 */
	public double getLearningTerm(){
		return (1-alpha)*(confidenceLevel*divergenceWithPreferences + (1-confidenceLevel)*divergenceWithUniformDistribution);
	}
	
	/**
	 * Get the current value of the Lagrangian term in objective function
	 * <br>(\lambda \times (\sum_{i=1}^{M}{w_i} -1) )
	 * @return
	 */
	public double getLagrangianTerm(){
		return lagrangian*(ArraysUtilities.getSumArray(learningVector) - 1);
	}
	
	/**
	 * Get the current value of the objective function
	 * @return
	 */
	public double getObjectiveFunction() {
		return m_CurrentObjective;
	}

	/**
	 * Set the vector of learned weights
	 * @param a A vector of weights
	 */
	public void setLearningVector(double[] a){
		this.learningVector = a.clone();
	}

	/**
	 * Get the vector of learned weights
	 * @return
	 */
	public double[] getLearningVector(){
		return this.learningVector;
	}
	
	/**
	 * Get the preferences vector
	 * @return
	 */
	public double[] getPreferencesVector(){
		return this.preferencesVector;
	}
	

	/**
	 * Get the number of iterations until convergence
	 * @return
	 */
	public int getNbIter(){return this.nbiter;}
	
	/**
	 * Get an array of clusters effectives
	 * @return
	 */
	public int[] getClustersEffective(){return this.clustersEffective.clone();}
	
	
	/**
	 * Get the data partition
	 * @return
	 */
	public int[] getPartition(){
		//trier les clusters par effectifs
		return this.partition.clone();
	}
	
	/**
	 * Get the number of discovered clusters
	 * @return
	 */
	public int getNbClusters(){
		List<Integer> map = new ArrayList<Integer>();
		for (int i = 0; i < this.partition.length; i++){
			if (!map.contains(partition[i])) map.add(partition[i]);
		}
		return map.size();
	}
	
	/**
	 * Get the computed centroids
	 * @return
	 */
	public List<double[]> getCentres() {
		// TODO Auto-generated method stub
		return this.centers;
	}
	
	
	/***
	 * add a point into cluster c
	 * @param centers cluster centers 
	 * @param point the data point
	 * @param c the id of the cluster
	 */
	private void add(List<double[]> centers, double[] point, int c){
		double[] tmp = centers.get(c);
		for (int i = 0; i < tmp.length; i++) tmp[i] = tmp[i] + point[i];
		centers.set(c,tmp);
	}	
	private void div(List<double[]> centers, int effectif, int c){
		double[] tmp = centers.get(c);
		for (int i = 0; i < tmp.length; i++)tmp[i] = tmp[i] / effectif;
		centers.set(c,tmp);
	}
	
	

	/* ***********  result traitement ************** */
	
	
	
	
	public int[] sortClustersByEffectifs(){
		int[] sorted_index = sortAndGetNewMappingIndices(clustersEffective);
		
		for(int i=0; i<partition.length; i++){
			partition[i] = sorted_index[partition[i]];
		}
		
		return partition;
	}
	
	/**
	 * Sort a table ascending the get the mapping indices of points in the new table
	 * @param table
	 * @return
	 */
	private int[] sortAndGetNewMappingIndices(int table[]) {
		int tableLength = table.length;
		int tempValue = 0;
		boolean permut;
		int[] newi = new int[tableLength];
		for (int j=0; j<tableLength; j++) newi[j] = j;
	
		do {
			
			permut = false;
			for (int i = 0; i < tableLength - 1; i++) {
				int[] index = newi.clone();
				
				if (table[i] > table[i + 1]) {
					// swap position
					tempValue = table[i];
					table[i] = table[i + 1];
					table[i + 1] = tempValue;
					newi[i]=index[i+1];
					newi[i+1]=index[i];
					permut = true;
				}
				
			}
		} while (permut);
		
		int[] mapping = new int[tableLength];
		for (int i=0; i<tableLength; i++){
			for (int j=0; j<tableLength; j++){
				if(i == newi[j]){
					mapping[i] = j;
				}
			}
		}
		return mapping;
	}
}
