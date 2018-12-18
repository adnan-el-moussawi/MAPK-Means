package comparator;

/**
 * Interface contains the definition of distance measures to compare data points 
 * @author Adnan EL MOUSSAWI
 *
 * @param <D> A data point is represented by a numeric value or an array of values 
 */
public interface Comparator<D> {
	/**
	 * Compute a distance between two points
	 * @param dat1 first data point
	 * @param dat2 second second data point
	 * @return The distance between two data points
	 */
	public double distance(D dat1, D dat2);
	
	/**
	 * Compute a squared distance
	 * @param dat1 first data point
	 * @param dat2 second second data point
	 * @return The distance between two data points
	 */
	public double distance2(D dat1, D dat2);
	
	
	/**
	 * Compute a weighted distance between two points
	 * @param dat1 first data point
	 * @param dat2 second second data point
	 * @param weigths generally a vector of weights of data attributes
	 * @return The distance between two data points
	 */
	public double weightedDistance(D dat1, D dat2, D weigths);
	
	/**
	 * Compute a squared weighted distance between two points
	 * @param dat1 first data point
	 * @param dat2 second second data point
	 * @param weigths generally a vector of weights of data attributes
	 * @return The distance between two data points
	 */
	public double weightedDistance2(D dat1, D dat2, D weights);
}
