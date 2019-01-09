package Main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import utilities.ArraysUtilities;

/**
 * This class contains several methods to compute an approximation of an optimal 
 * vector of weights on the attributes of the data, basing on some supervised measures 
 * on an exiting ground truth partition. (Adnan et al. 2016)<br>
 * 
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
 * 
 * @author Adnan EL MOUSSAWI
 *
 */
public class ComputeWeights {

	/**
	 * Compute the vector of optimal weights basing on the intra-cluster distances.
	 * <br>(the higher the distance on a dimension is, the higher is the weight of this dimension)
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return A vector of weights
	 */
	public static double[] computeOptimalWeightsIntra(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
		
		int nbDim = normalizedData.get(0).length;
		
		
		double[] weights = new double[nbDim],
				intra = computeIntra(normalizedData, groundTruthPartition, k),
				v = computeOptimlaWeightsVarianceBased(normalizedData, groundTruthPartition, k);

		for (int d=0; d<nbDim; d++){
			weights[d] =  v[d]/intra[d];
		}
		weights = ArraysUtilities.normalizeArray(weights);
		
		return weights;
	}


	/**
	 * Compute a vector of optimal weights basing on the inter-clusters distance.
	 * <br>(the higher the distance on a dimension is, the higher is the weight of this dimension)
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @param cp
	 * @return A vector of weights
	 */
	public static double[] computeOptimalWeightsInterClusterBased(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
				
		double[] weights = computeInter(normalizedData, groundTruthPartition, k);
		double sum = ArraysUtilities.getSumArray(weights);	
		for (int d=0; d<weights.length; d++){
			weights[d] =  weights[d]/sum;
		}
		
		return weights;
	}

	/**
	 * Compute a vector of optimal weights using the Fisher index of ANOVA-test.
	 * <br>(the higher the correlation between the dimension and the ground truth is, the higher is its weight)
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return A vector of weights
	 */
	public static double[] computeOptimalWeightAnova(List<double[]> normalizedData, List<String> groundTruthPartition,
			int k) {
		int nbDim = normalizedData.get(0).length,
				nbData = normalizedData.size();
		
		
		double[] weights = new double[nbDim],
				intra = computeIntra(normalizedData, groundTruthPartition, k),
				inter = computeInter(normalizedData, groundTruthPartition, k);
	
		for (int d=0; d<nbDim; d++){
			weights[d] =  ( (nbData-k)*inter[d])/( (k-1)*intra[d]);
		}
		weights = ArraysUtilities.normalizeArray(weights);
		
		return weights;
	}


	/**
	 * Compute a vector of irrelevant weights using the Fisher index of ANOVA-test.
	 * <br>(the less the correlation between the dimension and the ground truth is, the higher is its weight)
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return A vector of weights
	 */
	public static double[] computeBadWeightAnova(List<double[]> normalizedData, List<String> groundTruthPartition,
			int k) {
		int nbDim = normalizedData.get(0).length,
				nbData = normalizedData.size();
		
		
		double[] weights = new double[nbDim],
				intra = computeIntra(normalizedData, groundTruthPartition, k),
				inter = computeInter(normalizedData, groundTruthPartition, k);
	
		for (int d=0; d<nbDim; d++){
			weights[d] =  ( (k-1)*intra[d])/( (nbData-k)*inter[d]);
		}
		weights = ArraysUtilities.normalizeArray(weights);
		
		return weights;
	}



	/**
	 * Compute the vector of optimal weights using the inverse of cluster distortion, as explained in (Sun et al., 2010) 
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return A vector of weights
	 */
	public static double[] computeOptimalWeightsSun(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
		
		int nbData = normalizedData.size(),
				nbDim = normalizedData.get(0).length;
		
		/* compute global center */
		double[] globaleCenter = ArraysUtilities.getAverage(normalizedData, nbDim);
		
		double[] v = new double[nbDim],	distortion = new double[nbDim],	inverse_distortion = new double[nbDim],
				weights = new double[nbDim];
		
		Arrays.fill(v, 0.0);
		
		/* compute dimension variance */
		for(int d=0; d<nbDim; d++){
			for(int i=0; i<nbData; i++){
				v[d] += Math.pow((normalizedData.get(i)[d] - globaleCenter[d]), 2.0);
			}
		}
		
		
		/* get clusters */
		String[] clusters = new HashSet<String>(groundTruthPartition).toArray(new String[k]);
		
		/* Compute the centroids */
		ArrayList<double[]> centers = new ArrayList<double[]>();		
		int[] effectives = new int[k];
		for(int c=0; c<k; c++){
			double[] tmp = new double[nbDim];
			Arrays.fill(tmp, 0.0);
			for (int i = 0; i < nbData; i++){
				if(groundTruthPartition.get(i).equals(clusters[c])){
					for (int d=0; d<nbDim; d++){
						tmp[d] += normalizedData.get(i)[d];
					}
					effectives[c]++;
				}
			}
			for (int d=0; d<nbDim; d++) tmp[d] = tmp[d]/effectives[c];
			centers.add(tmp.clone());
		}
		
		/* compute intra-class distortion on each dimension */
		Arrays.fill(distortion, 0.0);
		double distortion_totale = 0.0;
		for (int d=0; d<nbDim; d++){
			for(int i=0; i<nbData; i++){
				for(int c=0; c<k; c++){
					//if(groundTruthPartition.get(c).equals(clusters[c])){
					if(groundTruthPartition.get(i).equals(clusters[c])){
						distortion[d] += Math.pow((normalizedData.get(i)[d] - centers.get(c)[d]), 2.0);
					}
				}
			}
			distortion[d] = distortion[d]/v[d];
			distortion_totale += distortion[d];
		}
	
		//System.out.println(toStringArray(distortion));
		
		/* compute inverse intra-class distortion on each dimension */
		Arrays.fill(inverse_distortion, 0.0);
		double inverse_distortion_totale=0.0;//to normalize such \sum weights = 1
		for (int d=0; d<nbDim; d++){
			inverse_distortion[d] = (distortion_totale)/distortion[d] -1.0;
			inverse_distortion_totale += inverse_distortion[d];
		}
		//System.out.println(toStringArray(inverse_distortion));
		
		/* compute normalized weights */
		Arrays.fill(weights, 0.0);
		for (int d=0; d<nbDim; d++){
			weights[d] =  inverse_distortion[d]/inverse_distortion_totale;
		}
		
		return weights;
	}


	/**
	 * Compute the vector of optimal weights using the inverse of inter-cluster distance 
	 * normalized by the variance on the dimension
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return A vector of weights
	 */
	public static double[] computeOptimlaWeightsVarianceBased(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
		
		int nbData = normalizedData.size(),
				nbDim = normalizedData.get(0).length;
		
		/* compute global center */
		double[] globaleCenter = ArraysUtilities.getAverage(normalizedData, nbDim);
		
		double[] weights = new double[nbDim];
		
		Arrays.fill(weights, 0.0);
		
		/* compute dimension variance */
		for(int d=0; d<nbDim; d++){
			for(int i=0; i<nbData; i++){
				weights[d] += Math.pow((normalizedData.get(i)[d] - globaleCenter[d]), 2.0);
			}
		}
		
		return weights;
	}
	



	/**
	 * Compute the inter-clusters distance
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return The vector of inter-clusters distance values on each dimension
	 */
	private static double[] computeInter(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
		
		int nbData = normalizedData.size(),
				nbDim = normalizedData.get(0).length;
		
		/* compute global center */
		double[] globaleCenter = ArraysUtilities.getAverage(normalizedData, nbDim);
		
		double[] inter = new double[nbDim];
		
		
		/* get clusters */
		String[] clusters = new HashSet<String>(groundTruthPartition).toArray(new String[k]);
		
		/* Compute the centroids */
		ArrayList<double[]> centers = new ArrayList<double[]>();
		int[] effectives = new int[k];
		for(int c=0; c<k; c++){
			double[] tmp = new double[nbDim];
			Arrays.fill(tmp, 0.0);
			for (int i = 0; i < nbData; i++){
				if(groundTruthPartition.get(i).equals(clusters[c])){
					for (int d=0; d<nbDim; d++){
						tmp[d] += normalizedData.get(i)[d];
					}
					effectives[c]++;
				}
			}
			for (int d=0; d<nbDim; d++) tmp[d] = tmp[d]/effectives[c];
			
			centers.add(tmp.clone());
		}
		
		Arrays.fill(inter, 0.0);
		/* compute inter class on each dimension */
		for (int d=0; d<nbDim; d++){
			for(int c=0; c<k; c++){
				inter[d] += effectives[c]*Math.pow((centers.get(c)[d] - globaleCenter[d]), 2.0);
			}
		}
		
		return inter.clone();
	}
	



	/**
	 * Compute the total intra-cluster distance
	 * @param normalizedData data set normalized between 0 and 1
	 * @param groundTruthPartition the ground truth partition of the data
	 * @param k number of cluster
	 * @return The vector of the total intra-cluster distance values on each dimension
	 */
	private static double[] computeIntra(List<double[]> normalizedData,
		List<String> groundTruthPartition, int k) {
		
		int nbData = normalizedData.size(),
				nbDim = normalizedData.get(0).length;
		
		
		double[] intra = new double[nbDim];
				
		/* get clusters */
		String[] clusters = new HashSet<String>(groundTruthPartition).toArray(new String[k]);
		
		/* Compute the centroids */
		ArrayList<double[]> centers = new ArrayList<double[]>();
		int[] effectives = new int[k];
		for(int c=0; c<k; c++){
			double[] tmp = new double[nbDim];
			Arrays.fill(tmp, 0.0);
			for (int i = 0; i < nbData; i++){
				if(groundTruthPartition.get(i).equals(clusters[c])){
					for (int d=0; d<nbDim; d++){
						tmp[d] += normalizedData.get(i)[d];
					}
					effectives[c]++;
				}
			}
			for (int d=0; d<nbDim; d++) tmp[d] = tmp[d]/effectives[c];
			
			centers.add(tmp.clone());
		}		
		
		
		/* compute intra-class on each dimension */		
		Arrays.fill(intra, 0.0);
		for (int d=0; d<nbDim; d++){
			for(int i=0; i<nbData; i++){
				for(int c=0; c<k; c++){
					//if(groundTruthPartition.get(c).equals(clusters[c])){
					if(groundTruthPartition.get(i).equals(clusters[c])){
						intra[d] += Math.pow((normalizedData.get(i)[d] - centers.get(c)[d]), 2.0);
					}
				}
			}
		}
		
		return intra.clone();
	}
	
	/**
	 * Adjust a given vector of weights such that the sum of weights >= 1/M 
	 * represents 90% of the total weight and 10% for the others
	 * @param weights A vector of weights
	 * @return An adjusted vector of weights
	 */
	public static double[] adjustWeights(double[] weights){
		int nbDim = weights.length;
		double sum_higher_weights = 0.0, sum_lower_weights = 0.0;
		double[] adjusted_weights = new double[nbDim];
		
		ArrayList<Integer> top3Indices = getTopIndices(weights);
		for(int i=0; i<nbDim; i++){
			if(top3Indices.contains(i)){
				sum_higher_weights += weights[i];
			}
			else{
				sum_lower_weights += weights[i];
			}
		}
		
		for(int i=0; i<nbDim; i++){
			if(top3Indices.contains(i)){
				adjusted_weights[i] = (weights[i]*0.9) / sum_higher_weights ;
			}
			else{
				adjusted_weights[i] = (weights[i]*0.1) / sum_lower_weights;
			}
		}
		
		
		return adjusted_weights.clone();
	}


	/**
	 * For a given vector of weights, get the top-n dimensions (by their indices) with the higher weights.
	 * <br>(n=1 if nbDim <= 3, n=2 if 3 < nbDim <= 6, n=3 if 6 < nbDim <= 12, else n=4)
	 * @param weights A vector of weights
	 * @return The vector of indices of the top-n dimensions
	 */
	private static ArrayList<Integer> getTopIndices(double weights[]) {
		double[] sorted_weights = weights.clone();
		int t_Length = weights.length;
		double tmp = 0.0;
		boolean permut;
		/* a vector the save the indices of the sorted dimensions */
		int[] sorted_Index = new int[t_Length];
		for (int j=0; j<t_Length; j++) sorted_Index[j] = j;
		
		/* sort the vector of weights */
		do {
			
			permut = false;
			for (int i = 0; i < t_Length - 1; i++) {
				int[] index = sorted_Index.clone();
				
				if (sorted_weights[i] > sorted_weights[i + 1]) {
					
					tmp = sorted_weights[i];
					sorted_weights[i] = sorted_weights[i + 1];
					sorted_weights[i + 1] = tmp;
					sorted_Index[i]=index[i+1];
					sorted_Index[i+1]=index[i];
					permut = true;
				}
				
			}
		} while (permut);
		
		ArrayList<Integer> l = new ArrayList<>();
		l.add(sorted_Index[t_Length-1]);
		
		if(t_Length > 3){
			l.add(sorted_Index[t_Length-2]);
		}
		if(t_Length > 6){
			l.add(sorted_Index[t_Length-3]);
		}
		
		if(t_Length > 12){
			l.add(sorted_Index[t_Length-4]);
		}
		return l;
	}
	
}
