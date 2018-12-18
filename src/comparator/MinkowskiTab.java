package comparator;

/**
 * This class implements the Minkowski distance
 * @author Adnan EL MOUSSAWI
 *
 */
public class MinkowskiTab implements Comparator<double[]> {
	/**
	 * The order of the Minkowski distance
	 */
	int order;
	
	/**
	 * @param order : if order = 1 than Manhattan distance, if order = 2 than euclidean distance
	 */
	public MinkowskiTab(int order){
		this.order = order;
	}
	
	/**
	 * Compute the distance between two data points
	 * @param v1 first data point
	 * @param v2 second data point
	 * @return distance between two vectors of double
	 */
	public double distance(double[] v1, double[] v2) {
		
		double somme = 0.0;
		for (int i = 0; i < v1.length; i++){
			somme += Math.pow(v1[i] - v2[i], order);
		}
		return Math.pow(somme, 1.0/order);
	}
	

	/**
	 * Compute the squared distance between two points
	 * @param v1 first data point
	 * @param v2 second data point
	 * @return squared distance between two vectors of double
	 */
	public double distance2(double[] v1, double[] v2) {
		return  Math.pow(distance(v1, v2), 2);
	}
	
	
	/**
	 * Compute the weighted distance between two points
	 * @param v1 first data point
	 * @param v2 second data point
	 * @param a vector of weights
	 * @return weighted distance between two vectors of double
	 */
	public double weightedDistance(double[] v1, double[] v2, double[] a) {
		double somme = 0.0;
		for (int i = 0; i < v1.length; i++){
			//System.out.println(a[i]*(v1[i] - v2[i]));
			somme += a[i]*Math.pow((v1[i] - v2[i]) , order);
		}
		//System.out.println(somme);
		return Math.pow(somme, 1.0/order);
	}
	
	
	/**
	 * Compute the squared weighted distance between two points
	 * @param v1 first data point
	 * @param v2 second data point
	 * @param a vector of weights
	 * @return squared weighted distance between two vectors of double
	 */
	public double weightedDistance2(double[] v1, double[] v2, double[] a) {
		return Math.pow(weightedDistance(v1, v2, a), 2);
	}
	
}
