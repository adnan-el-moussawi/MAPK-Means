package evaluation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Some External indices for clustering evaluation
 * @author Adnan EL MOUSSAWI
 *
 */

public class ExternalEvaluation{
	
	/**
	 * Compute the Adjusted Rand Index de (Hubert, L. & Arabie, P. Comparing partitions Journal of Classification, 1985, 2, 193-218)
	 * <br>(The compared partitions must have the same size)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as an int array)
	 * @return The ARI value between two partitions
	 */
	public static float getARI(int[] partition1, int[] partition2) {
		// TODO Auto-generated method stub
		List<String> p2 = new ArrayList<String>();
		for (int value : partition2){
			p2.add(value+"");
		}
		return getARI(partition1, p2);
	}
	
	/**
	 * Compute the Adjusted Rand Index de (Hubert, L. & Arabie, P. Comparing partitions Journal of Classification, 1985, 2, 193-218)
	 * <br>(The compared partitions must have the same size)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return The ARI value between two partitions
	 */
	public static float getARI(int[] partition1, List<String> partition2){
		int[][] contingency = contingency(partition1, partition2);
		
		float nxy = 0, nxy2 = 0;
		float[] a = new float[contingency.length], b = new float[contingency[0].length];
		for(int i=0; i<a.length; i++){
			for(int j=0; j<b.length; j++){
				a[i] += contingency[i][j];
				b[j] += contingency[i][j];
				nxy += contingency[i][j];
				nxy2 += contingency[i][j]*contingency[i][j];
			}
		}

		float nx = 0, nx2 = 0;
		for (int i=0; i < a.length; i++){
			nx += a[i];
			nx2 += a[i]*a[i];
		}

		float ny = 0, ny2 = 0;
		for (int j=0; j < b.length; j++){
			ny += b[j];
			ny2 += b[j]*b[j];
		}


		// Expected index
		float ARI = 0;
		float n = 0,sumxy=0,sumx=0,sumy=0;

		n = ny;

		sumxy = nxy2 - nxy;
		sumx = nx2 - nx;
		sumy = ny2 - ny;

		float mulxy = sumx * sumy / (n*(n-1)); 

		ARI = (sumxy - mulxy)/( (sumx+sumy)/2 - mulxy);
		return ARI;
	}

	
	/**
	 * 
	 * The approximate ARI formula as used by Milligan and introduced in Morrey and Agresti in 1984
	 * <br>(The compared partitions must have the same size)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as an int array)
	 * @return The approximate ARI value between two partitions
	 */
	public static float getApproximateARI(int[] partition1, int[] partition2) {
		// TODO Auto-generated method stub
		List<String> p2 = new ArrayList<String>();
		for (int value : partition2){
			p2.add(value+"");
		}
		return getApproximateARI(partition1, p2);
	}



	/**
	 * 
	 * The approximate ARI formula as used by Milligan and introduced in Morrey and Agresti in 1984
	 * <br>(The compared partitions must have the same size)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return The approximate ARI value between two partitions
	 */
	public static float getApproximateARI(int[] partition1, List<String> partition2){
		int[][] contingency = contingency(partition1, partition2);
		
		float nxy2 = 0;
		float[] a = new float[contingency.length], b = new float[contingency[0].length];
		for(int i=0; i<a.length; i++){
			for(int j=0; j<b.length; j++){
				a[i] += contingency[i][j];
				b[j] += contingency[i][j];
				nxy2 += contingency[i][j]*contingency[i][j];
			}
		}

		float nx2 = 0;
		for (int i=0; i < a.length; i++){
			nx2 += a[i]*a[i];
		}

		float ny = 0, ny2 = 0;
		for (int j=0; j < b.length; j++){
			ny += b[j];
			ny2 += b[j]*b[j];
		}
		
		float ARI = 0;
		float n = 0,sumxy=0,sumx=0,sumy=0;
		
		n = ny;
		sumxy = nxy2;
		sumx = nx2;
		sumy = ny2;
		
		float mulxy = sumx * sumy / (n*n); 
		
		ARI = (sumxy - mulxy)/( (sumx+sumy)/2 - mulxy);
		return ARI;
	}
	

	/**
	 * Compute the contingency matrix between two partitions
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as an int array)
	 * @return Contingency matrix
	 */
	static protected int[][] contingency(int[] partition1, int[] partition2){
		List<String> p2 = new ArrayList<String>();
		for (int value : partition2){
			p2.add(value+"");
		}
		return contingency(partition1, p2);
	}



	/**
	 * Compute the contingency matrix between two partitions
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return Contingency matrix
	 */
	static protected int[][] contingency(int[] partition,  List<String> partition2){
		//int nbcluster1 = Evaluation.getNbClusters(partition);
		//int nbcluster2 = Evaluation.getNbClusters2(partition2);
		
		// Construction a dictionary of label for each partition :
		// to process empty clusters and convert string labels into int value
		Map<Integer, Integer> dico1 = new HashMap<Integer, Integer>();
		Map<String, Integer> dico2 = new HashMap<String, Integer>();
		int count = 0;
		for (int value : partition){
			if (!dico1.containsKey(value)){
				dico1.put(value, count);
				count ++;
			}
		}
		count = 0;
		for (String value : partition2){
			if (!dico2.containsKey(value)){
				dico2.put(value, count);
				count ++;
			}
		}
		int nbcluster1 = dico1.size();
		int nbcluster2 = dico2.size();
		
		int[][] contingency = new int [nbcluster1][nbcluster2];
		
		for (int i = 0; i < partition.length; i++){
			contingency[dico1.get(partition[i])][dico2.get(partition2.get(i))]++;
		}
		
		return contingency;
	}



	/**
	 * Compute the contingency matrix between two partitions
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return Contingency matrix
	 */
	static protected int[][] contingency(List<String> partition1,  List<String> partition2){
		//int nbcluster1 = Evaluation.getNbClusters(partition);
		//int nbcluster2 = Evaluation.getNbClusters2(partition2);
		
		// Construction a dictionary of label for each partition :
		// to process empty clusters and convert string labels into int value
		Map<String, Integer> dico1 = new HashMap<String, Integer>();
		Map<String, Integer> dico2 = new HashMap<String, Integer>();
		int count = 0;
		for (String value : partition1){
			if (!dico1.containsKey(value)){
				dico1.put(value, count);
				count ++;
			}
		}
		count = 0;
		for (String value : partition2){
			if (!dico2.containsKey(value)){
				dico2.put(value, count);
				count ++;
			}
		}
		int nbcluster1 = dico1.size();
		int nbcluster2 = dico2.size();
		
		int[][] contingency = new int [nbcluster1][nbcluster2];
		
		for (int i = 0; i < partition1.size(); i++){
			contingency[dico1.get(partition1.get(i))][dico2.get(partition2.get(i))]++;
		}
		
		return contingency;
	}



	/**
	 * Compute the joint-normalized mutual information (Sun et al., 2010, Clustering with Feature Order Preferences,
	 * Intelligent Data Analysis, vol. 14, no. 4, pp. 479-495, 2010)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as an int array)
	 * @return Joint Normalized Mutual information
	 */
	public static double entropyNMISun2010(int[] partition1, int[] partition2) {
		int[][] contingency = contingency(partition1, partition2);
		double n = partition1.length;

		double[] a = new double[contingency.length], b = new double[contingency[0].length];
		Arrays.fill(a, 0.0); Arrays.fill(b, 0.0);
		
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				a[i] += contingency[i][j];
				b[j] += contingency[i][j];
			}
		}
		
		double entropyFirst = 0.0;
		// iterate over first clustering
		for (int i = 0; i < a.length; i++) {
			if (a[i] > 0) {
				entropyFirst += a[i] * Math.log(a[i] / n);
			}
		}

		double entropySecond = 0.0;
		// iterate over first clustering
		for (int j = 0; j < b.length; j++) {
			if (b[j] > 0) {
				entropySecond += b[j] * Math.log(b[j] / n);
			}
		}
		
		double entropyJoint = 0.0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				if(contingency[i][j]>0){
					entropyJoint += contingency[i][j] * Math.log( (n*contingency[i][j])/(a[i]*b[j]) );
				}
			}
		}
		
		if (entropyJoint == 0.0) {
			return 0.0;
		}
		return ( entropyJoint / Math.sqrt(entropyFirst*entropySecond) );
	}



	/**
	 * Compute the joint-normalized mutual information (Sun et al., 2010, Clustering with Feature Order Preferences,
	 * Intelligent Data Analysis, vol. 14, no. 4, pp. 479-495, 2010)
	 * @param partition1 Discovered partition (given as a list of String)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return Joint Normalized Mutual information
	 */
	public static double entropyNMISun2010(List<String> partition1, List<String> partition2) {
		int[][] contingency = contingency(partition1, partition2);
		double n = partition1.size();
	
		double[] a = new double[contingency.length], b = new double[contingency[0].length];
		Arrays.fill(a, 0.0); Arrays.fill(b, 0.0);
		
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				a[i] += contingency[i][j];
				b[j] += contingency[i][j];
			}
		}
		
		double entropyFirst = 0.0;
		// iterate over first clustering
		for (int i = 0; i < a.length; i++) {
			if (a[i] > 0) {
				entropyFirst += a[i] * Math.log(a[i] / n);
			}
		}
	
		double entropySecond = 0.0;
		// iterate over first clustering
		for (int j = 0; j < b.length; j++) {
			if (b[j] > 0) {
				entropySecond += b[j] * Math.log(b[j] / n);
			}
		}
		
		double entropyJoint = 0.0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				if(contingency[i][j]>0){
					entropyJoint += contingency[i][j] * Math.log( (n*contingency[i][j])/(a[i]*b[j]) );
				}
			}
		}
		
		if (entropyJoint == 0.0) {
			return 0.0;
		}
		return ( entropyJoint / Math.sqrt(entropyFirst*entropySecond) );
	}

	/**
	 * Compute the joint-normalized mutual information (Sun et al., 2010, Clustering with Feature Order Preferences,
	 * Intelligent Data Analysis, vol. 14, no. 4, pp. 479-495, 2010)
	 * @param partition1 Discovered partition (given as an int array)
	 * @param partition2 Ground truth partition (given as a list of String)
	 * @return Joint Normalized Mutual information
	 */
	public static double entropyNMISun2010(int[] partition1, List<String> partition2) {
		int[][] contingency = contingency(partition1, partition2);
		double n = partition1.length;
	
		double[] a = new double[contingency.length], b = new double[contingency[0].length];
		Arrays.fill(a, 0.0); Arrays.fill(b, 0.0);
		
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				a[i] += contingency[i][j];
				b[j] += contingency[i][j];
			}
		}
		
		double entropyFirst = 0.0;
		// iterate over first clustering
		for (int i = 0; i < a.length; i++) {
			if (a[i] > 0) {
				entropyFirst += a[i] * Math.log(a[i] / n);
			}
		}
	
		double entropySecond = 0.0;
		// iterate over first clustering
		for (int j = 0; j < b.length; j++) {
			if (b[j] > 0) {
				entropySecond += b[j] * Math.log(b[j] / n);
			}
		}
		
		double entropyJoint = 0.0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				if(contingency[i][j]>0){
					entropyJoint += contingency[i][j] * Math.log( (n*contingency[i][j])/(a[i]*b[j]) );
				}
			}
		}
		
		if (entropyJoint == 0.0) {
			return 0.0;
		}
		return ( entropyJoint / Math.sqrt(entropyFirst*entropySecond) );
	}
}
