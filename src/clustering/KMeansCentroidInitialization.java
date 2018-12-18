package clustering;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import comparator.Comparator;
import comparator.MinkowskiTab;
import utilities.ArraysUtilities;

/**
 * @author Adnan EL MOUSSAWI
 *
 */
public class KMeansCentroidInitialization {


	private Random r = new Random();
	private Comparator<double[]> cp = new MinkowskiTab(2);

	/***
	 * This method chooses k centroids random from data 
	 * @param data : the data set ({@link List[<code>double</code>]})
	 * @param k : number of clusters
	 * @return a list of k centroids
	 */
	public List<double[]> initializationRandom(List<double[]> data, int k){
		List<Integer> indices = new ArrayList<Integer>();
		int candidat = -1;
		for (int i = 0; i < k; i++){
			do {
				candidat = r.nextInt(data.size());
			} while (indices.contains(candidat));
			indices.add(candidat);
		}
		List<double[]> centres = new ArrayList<double[]>();
		for (int i = 0; i < k; i++){
			centres.add(data.get(indices.get(i)));
		}
		return centres;
	}
	
	/***
	 * initialization de k centres avec les k points de donnees les plus eloignes et C0 = X0
	 * @param cp {@link comparator.Comparator<double[]>}
	 * @param data : the data set ({@link List[<code>double</code>]})
	 * @param k : number of clusters
	 * @return a list of k centroids
	 */
	public List<double[]> initializationKFarthestFix(List<double[]> data, int k){
		List<double[]> centres = new ArrayList<double[]>();
		List<double[]> tempData = new ArrayList<double[]>(data);
		
		centres.add(tempData.get(0));
		tempData.remove(0);
		
		double[] pt; double minDc, tmp, maxD; int index = -1;
		
		while (centres.size() < k){
			maxD = -1.0;
			
			for (int i=0; i<tempData.size(); i++){
				pt = tempData.get(i).clone();
				minDc = cp.distance2(pt, centres.get(0));
				for (int c=1; c<centres.size(); c++){
					tmp = cp.distance2(pt, centres.get(c));
					if (tmp < minDc)
						minDc = tmp;
				}
				if (i == 0){
					maxD = minDc;
					index = 0;
				}
				else {
					if (maxD < minDc){
						maxD = minDc;
						index = i;
					}
				}
			}
			
			centres.add(tempData.get(index));
			tempData.remove(index);
		}
		return centres;
	}
	
	/***
	 * initialization de k centres avec les k points de donnees les plus eloignes et C0 = X0
	 * @param cp {@link comparator.Comparator<double[]>}
	 * @param data : the data set ({@link List[<code>double</code>]})
	 * @param k : number of clusters
	 * @return a list of k centroids
	 */
	public List<double[]> initializationKFarthestFixOptimized(List<double[]> data, int k){
		List<double[]> centres = new ArrayList<double[]>();
		List<double[]> tempData = new ArrayList<double[]>(data);
		
		int nbDimension = data.get(0).length;
	   
		/* compute global center */
		double[] globaleCenter = ArraysUtilities.getAverage(data, nbDimension);
		
		double average_gdistance = 0.0,
				sd_globale_distance = 0.0, inf_gdistance, sup_gdistance;
		double[] distances = new double[data.size()];
		/* compute distance average */
		for(int i=0; i<data.size(); i++){
			distances[i] = cp.distance2(data.get(i), globaleCenter);
			//System.out.println(distances[i]);
			average_gdistance += distances[i];
		}
		average_gdistance = average_gdistance/data.size();
		sd_globale_distance = ArraysUtilities.getSdArray(distances);

		inf_gdistance = average_gdistance - 1*sd_globale_distance;
		sup_gdistance = average_gdistance + 1*sd_globale_distance;
		int tmpDataSize = tempData.size();
		for (int i=0; i<tmpDataSize; i++){
			double d_to_gc = cp.distance2(data.get(i), globaleCenter);
			if(d_to_gc < inf_gdistance || d_to_gc > sup_gdistance){
				tempData.remove(i);
				tmpDataSize--;
			}
		}
		
		centres.add(tempData.get(0));
		tempData.remove(0);
		
		double[] pt; double minDc, tmp, maxD; int index = -1;
		
		while (centres.size() < k){
			maxD = -1.0;
			
			for (int i=0; i<tempData.size(); i++){
				pt = tempData.get(i).clone();
				minDc = cp.distance2(pt, centres.get(0));
				for (int c=1; c<centres.size(); c++){
					tmp = cp.distance2(pt, centres.get(c));
					if (tmp < minDc)
						minDc = tmp;
				}
				
				if (maxD < minDc){
					maxD = minDc;
					index = i;
				}
			}
			
			centres.add(tempData.get(index));
			tempData.remove(index);
		}
		return centres;
	}

	
	/***
	 * initialization de k centres avec les k points de donnees les plus eloignes et C0 est aleatoire
	 * @param cp {@link comparator.Comparator<double[]>}
	 * @param data : the data set ({@link List[<code>double</code>]})
	 * @param k : number of clusters
	 * @return a list of k centroids
	 */
	public List<double[]> initializationKFarthestRandom(List<double[]> data, int k){
		List<double[]> centres = new ArrayList<double[]>();
		List<double[]> tempData = new ArrayList<double[]>(data);
		int randomNumber = r.nextInt(data.size());
		
		centres.add(tempData.get(randomNumber));
		tempData.remove(randomNumber);
		
		double[] pt; double minDc, tmp, maxD; int index = -1;
		
		while (centres.size() < k){
			maxD = -1.0;
			
			for (int i=0; i<tempData.size(); i++){
				pt = tempData.get(i).clone();
				minDc = cp.distance2(pt, centres.get(0));
				for (int c=1; c<centres.size(); c++){
					tmp = cp.distance2(pt, centres.get(c));
					if (tmp < minDc)
						minDc = tmp;
				}
				if (i == 0){
					maxD = minDc;
					index = 0;
				}
				else {
					if (maxD < minDc){
						maxD = minDc;
						index = i;
					}
				}
			}
			
			centres.add(tempData.get(index));
			tempData.remove(index);
		}
		return centres;
	}
	
	/**

	 * Use K-means++ to choose the initial centers.
	 *
	 * @param <T> type of the points to cluster

	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers

	 */
	public List<double[]> initializationKMeansPlusPlus(List<double[]> data, final int k) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<double[]> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	
	    // Set the corresponding element in this array to indicate when
	    // elements of pointList are no longer available.
	    final boolean[] taken = new boolean[numPoints];
	
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
	
	    // Choose one center uniformly at random from among the data points.
	    // if center is outlier, choose another
	    int nbDimension = data.get(0).length;
	    double[] avg = ArraysUtilities.getAverage(data, nbDimension);
	    double[] sd = ArraysUtilities.getSd(data, nbDimension);
	    int tmp = 0;
	    boolean bad = true;
	    while(bad){
	    	bad=false;
	    	tmp = r.nextInt(numPoints);
	    	double[] tmpPoint = pointList.get(tmp);
	    	for(int j=0; j<nbDimension; j++){
	    		if( tmpPoint[j]>=(avg[j]+2*sd[j]) || tmpPoint[j]<=(avg[j]-2*sd[j])){
	    			bad = true;
	    		}
	    	}
	    }
	    
	    
	    final int firstPointIndex = tmp; 
	    final double[] firstPoint = pointList.get(firstPointIndex);
	
	    resultSet.add(firstPoint.clone());
	
	    // Must mark it as taken
	    taken[firstPointIndex] = true;
	
	    // To keep track of the minimum distance squared of elements of
	    // pointList to elements of resultSet.
	    final double[] minDistSquared = new double[numPoints];
	
	    // Initialize the elements.  Since the only point in resultSet is firstPoint,
	    // this is very easy.
	    for (int i = 0; i < numPoints; i++) {
	        if (i != firstPointIndex) { // That point isn't considered
	            double d = cp.distance2(firstPoint, pointList.get(i));
	            minDistSquared[i] = d*d;
	        }
	    }
	
	    while (resultSet.size() < k) {
	
	        // Sum up the squared distances for the points in pointList not
	        // already taken.
	        double distSqSum = 0.0;
	
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                distSqSum += minDistSquared[i];
	            }
	        }
	
	        // Add one new data point as a center. Each point x is chosen with
	        // probability proportional to D(x)2
	        final double rn = r.nextDouble() * distSqSum;
	
	        // The index of the next point to be added to the resultSet.
	        int nextPointIndex = -1;
	
	        // Sum through the squared min distances again, stopping when
	        // sum >= r.
	        double sum = 0.0;
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                sum += minDistSquared[i];
	                if (sum >= rn) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // If it's not set to >= 0, the point wasn't found in the previous
	        // for loop, probably because distances are extremely small.  Just pick
	        // the last available point.
	        if (nextPointIndex == -1) {
	            for (int i = numPoints - 1; i >= 0; i--) {
	                if (!taken[i]) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // We found one.
	        if (nextPointIndex >= 0) {
	
	            double[] p = pointList.get(nextPointIndex);
	
	            resultSet.add(p.clone());
	
	            // Mark it as taken.
	            taken[nextPointIndex] = true;
	
	            if (resultSet.size() < k) {
	                // Now update elements of minDistSquared.  We only have to compute
	                // the distance to the new center to do this.
	                for (int j = 0; j < numPoints; j++) {
	                    // Only have to worry about the points still not taken.
	                    if (!taken[j]) {
	                        double d = cp.distance2(p, pointList.get(j));
	                        double d2 = d * d;
	                        if (d2 < minDistSquared[j]) {
	                            minDistSquared[j] = d2;
	                        }
	                    }
	                }
	            }
	
	        } else {
	            // None found --
	            // Break from the while loop to prevent
	            // an infinite loop.
	            break;
	        }
	    }
	
	    return resultSet;
	}
	
	/**

	 * Use K-means++ to choose the initial centers, first center = first data point
	 *
	 * @param <T> type of the points to cluster

	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers

	 */
	protected List<double[]> initializationKMeansPlusPlusFix(List<double[]> data, final int k) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<double[]> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	
	    // Set the corresponding element in this array to indicate when
	    // elements of pointList are no longer available.
	    final boolean[] taken = new boolean[numPoints];
	
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
	    
	    
	    final int firstPointIndex = 0; 
	    final double[] firstPoint = pointList.get(firstPointIndex);
	
	    resultSet.add(firstPoint.clone());
	
	    // Must mark it as taken
	    taken[firstPointIndex] = true;
	
	    // To keep track of the minimum distance squared of elements of
	    // pointList to elements of resultSet.
	    final double[] minDistSquared = new double[numPoints];
	
	    // Initialize the elements.  Since the only point in resultSet is firstPoint,
	    // this is very easy.
	    for (int i = 1; i < numPoints; i++) {
	        double d = cp.distance2(firstPoint, pointList.get(i));
	        minDistSquared[i] = d*d;
	    }
	
	    while (resultSet.size() < k) {
	
	        // Sum up the squared distances for the points in pointList not
	        // already taken.
	        double distSqSum = 0.0;
	
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                distSqSum += minDistSquared[i];
	            }
	        }
	
	        // Add one new data point as a center. Each point x is chosen with
	        // probability proportional to D(x)2
	        final double rn = r.nextDouble() * distSqSum;
	
	        // The index of the next point to be added to the resultSet.
	        int nextPointIndex = -1;
	
	        // Sum through the squared min distances again, stopping when
	        // sum >= r.
	        double sum = 0.0;
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                sum += minDistSquared[i];
	                if (sum >= rn) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // If it's not set to >= 0, the point wasn't found in the previous
	        // for loop, probably because distances are extremely small.  Just pick
	        // the last available point.
	        if (nextPointIndex == -1) {
	            for (int i = numPoints - 1; i >= 0; i--) {
	                if (!taken[i]) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // We found one.
	        if (nextPointIndex >= 0) {
	
	            double[] p = pointList.get(nextPointIndex);
	
	            resultSet.add(p.clone());
	
	            // Mark it as taken.
	            taken[nextPointIndex] = true;
	
	            if (resultSet.size() < k) {
	                // Now update elements of minDistSquared.  We only have to compute
	                // the distance to the new center to do this.
	                for (int j = 0; j < numPoints; j++) {
	                    // Only have to worry about the points still not taken.
	                    if (!taken[j]) {
	                        double d = cp.distance2(p, pointList.get(j));
	                        double d2 = d * d;
	                        if (d2 < minDistSquared[j]) {
	                            minDistSquared[j] = d2;
	                        }
	                    }
	                }
	            }
	
	        } else {
	            // None found --
	            // Break from the while loop to prevent
	            // an infinite loop.
	            break;
	        }
	    }
	
	    return resultSet;
	}
	
	
	

	
	/***
	 * initialization de k centres avec les k premiers points
	 * @param data : the data set ({@link List[<code>double</code>]})
	 * @param k : number of clusters
	 * @return a list of k centroids
	 */
	public List<double[]> initializationFix(List<double[]> data, int k){
		List<double[]> centres = new ArrayList<double[]>();
		for (int i = 0; i < k; i++){
			centres.add(data.get(i));
		}
		return centres;
	}

	/**
	
	 * Use K-means++ to choose the initial centers.
	 *
	 * @param <T> type of the points to cluster
	
	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @return the initial centers
	
	 */
	protected List<double[]> initializationKMeansPlusPlusRecursive(List<double[]> data, final int k, int nbExecutions) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    //final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<double[]> (data));
	    List<double[]> pointList = new ArrayList<double[]> (data);
	
	    // The number of points in the list.
	    int numPoints = pointList.size();
	
	    // Set the corresponding element in this array to indicate when
	    // elements of pointList are no longer available.
	    boolean[] taken = new boolean[numPoints];
	
	    // The resulting list of initial centers.
	    List<double[]> resultSet = new ArrayList<double[]>();
	
	    // Choose one center uniformly at random from among the data points.
	    // if center is outlier, choose another
	    int nbDimension = data.get(0).length;
	    double[] avg = ArraysUtilities.getAverage(data, nbDimension);
	    double[] sd = ArraysUtilities.getSd(data, nbDimension);
	    int tmp = 0;
	    boolean bad = true;
	    while(bad){
	    	bad=false;
	    	tmp = r.nextInt(numPoints);
	    	double[] tmpPoint = pointList.get(tmp);
	    	for(int j=0; j<nbDimension; j++){
	    		if( tmpPoint[j]>=(avg[j]+2*sd[j]) || tmpPoint[j]<=(avg[j]-2*sd[j])){
	    			bad = true;
	    		}
	    	}
	    }
	    
	    
	    final int firstPointIndex = tmp; 
	    final double[] firstPoint = pointList.get(firstPointIndex);
	    
	    /*
	     * remove data points from pointList or data at each add
	     * use pointList or data as entry to recursive function
	     * with  nbExecutions--
	     */
	
	    resultSet.add(firstPoint.clone());
        data.remove(firstPoint.clone());
	
	    // Must mark it as taken
	    taken[firstPointIndex] = true;
	
	    // To keep track of the minimum distance squared of elements of
	    // pointList to elements of resultSet.
	    final double[] minDistSquared = new double[numPoints];
	
	    // Initialize the elements.  Since the only point in resultSet is firstPoint,
	    // this is very easy.
	    for (int i = 0; i < numPoints; i++) {
	        if (i != firstPointIndex) { // That point isn't considered
	            double d = cp.distance2(firstPoint, pointList.get(i));
	            minDistSquared[i] = d*d;
	        }
	    }
	
	    while (resultSet.size() < k) {
	
	        // Sum up the squared distances for the points in pointList not
	        // already taken.
	        double distSqSum = 0.0;
	
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                distSqSum += minDistSquared[i];
	            }
	        }
	
	        // Add one new data point as a center. Each point x is chosen with
	        // probability proportional to D(x)2
	        final double rn = r.nextDouble() * distSqSum;
	
	        // The index of the next point to be added to the resultSet.
	        int nextPointIndex = -1;
	
	        // Sum through the squared min distances again, stopping when
	        // sum >= r.
	        double sum = 0.0;
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                sum += minDistSquared[i];
	                if (sum >= rn) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // If it's not set to >= 0, the point wasn't found in the previous
	        // for loop, probably because distances are extremely small.  Just pick
	        // the last available point.
	        if (nextPointIndex == -1) {
	            for (int i = numPoints - 1; i >= 0; i--) {
	                if (!taken[i]) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // We found one.
	        if (nextPointIndex >= 0) {
	
	            double[] p = pointList.get(nextPointIndex);
	
	            resultSet.add(p.clone());
	            data.remove(p.clone());
	
	            // Mark it as taken.
	            taken[nextPointIndex] = true;
	
	            if (resultSet.size() < k) {
	                // Now update elements of minDistSquared.  We only have to compute
	                // the distance to the new center to do this.
	                for (int j = 0; j < numPoints; j++) {
	                    // Only have to worry about the points still not taken.
	                    if (!taken[j]) {
	                        double d = cp.distance2(p, pointList.get(j));
	                        double d2 = d * d;
	                        if (d2 < minDistSquared[j]) {
	                            minDistSquared[j] = d2;
	                        }
	                    }
	                }
	            }
	
	        } else {
	            // None found --
	            // Break from the while loop to prevent
	            // an infinite loop.
	            break;
	        }
	    }
	    if( (--nbExecutions)>0){
	    	resultSet.addAll(initializationKMeansPlusPlusRecursive(data, k, nbExecutions));
	    }
	    return resultSet;
	}

	/**
	
	 * Use K-means++ to choose the initial centers, first center = first data point
	 *
	 * @param <T> type of the points to cluster
	
	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers
	
	 */
	protected List<double[]> initializationKMeansPlusPlusFixRecursive(List<double[]> data, final int k, int nbExecutions) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<double[]> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	
	    // Set the corresponding element in this array to indicate when
	    // elements of pointList are no longer available.
	    final boolean[] taken = new boolean[numPoints];
	
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
	    
	    
	    final int firstPointIndex = 0; 
	    final double[] firstPoint = pointList.get(firstPointIndex);
	
	    resultSet.add(firstPoint.clone());
	
	    // Must mark it as taken
	    taken[firstPointIndex] = true;
	
	    // To keep track of the minimum distance squared of elements of
	    // pointList to elements of resultSet.
	    final double[] minDistSquared = new double[numPoints];
	
	    // Initialize the elements.  Since the only point in resultSet is firstPoint,
	    // this is very easy.
	    for (int i = 1; i < numPoints; i++) {
	        double d = cp.distance2(firstPoint, pointList.get(i));
	        minDistSquared[i] = d*d;
	    }
	
	    while (resultSet.size() < k) {
	
	        // Sum up the squared distances for the points in pointList not
	        // already taken.
	        double distSqSum = 0.0;
	
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                distSqSum += minDistSquared[i];
	            }
	        }
	
	        // Add one new data point as a center. Each point x is chosen with
	        // probability proportional to D(x)2
	        final double rn = r.nextDouble() * distSqSum;
	
	        // The index of the next point to be added to the resultSet.
	        int nextPointIndex = -1;
	
	        // Sum through the squared min distances again, stopping when
	        // sum >= r.
	        double sum = 0.0;
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                sum += minDistSquared[i];
	                if (sum >= rn) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // If it's not set to >= 0, the point wasn't found in the previous
	        // for loop, probably because distances are extremely small.  Just pick
	        // the last available point.
	        if (nextPointIndex == -1) {
	            for (int i = numPoints - 1; i >= 0; i--) {
	                if (!taken[i]) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // We found one.
	        if (nextPointIndex >= 0) {
	
	            double[] p = pointList.get(nextPointIndex);
	
	            resultSet.add(p.clone());
	
	            // Mark it as taken.
	            taken[nextPointIndex] = true;
	
	            if (resultSet.size() < k) {
	                // Now update elements of minDistSquared.  We only have to compute
	                // the distance to the new center to do this.
	                for (int j = 0; j < numPoints; j++) {
	                    // Only have to worry about the points still not taken.
	                    if (!taken[j]) {
	                        double d = cp.distance2(p, pointList.get(j));
	                        double d2 = d * d;
	                        if (d2 < minDistSquared[j]) {
	                            minDistSquared[j] = d2;
	                        }
	                    }
	                }
	            }
	
	        } else {
	            // None found --
	            // Break from the while loop to prevent
	            // an infinite loop.
	            break;
	        }
	    }
	
	    return resultSet;
	}

	/**
	 * Variance Partitioning Initialization,  
	 *
	 * @param <T> type of the points to cluster
	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers
	 */
	public List<double[]> VariancePartitioningInitialization(List<double[]> data, final int k) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	    int[] partitions = new int[numPoints];
	    Arrays.fill(partitions, 0);
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
		List<double[]> c1 = new ArrayList<double[]>(),
				c2 = new ArrayList<double[]>();
	    List < List<double[]> > clusters = new ArrayList< List<double[]> >();
	
	    // Choose one center uniformly at random from among the data points.
	    // if center is outlier, choose another
	    int nbDimension = data.get(0).length;
	    List<double[]> trainingData = new ArrayList<double[]>();
	    clusters.add(data);
	    
	    int bestCluster = 0, iVMax;
		double bestInter = -1.0;
		double intra, avg, maxVar;
		
		for(int c=1; c<k; c++){
	    	trainingData = clusters.get(bestCluster);
	    	
	    	clusters.remove(bestCluster);
	    	//System.out.println(trainingData.size());
	    	double[] center = ArraysUtilities.getAverage(trainingData, nbDimension);
	    	double[] intraVariance = computeIntraVariance(trainingData);
	    	maxVar = intraVariance[0];
	    	iVMax=0;
	    	for(int i=1; i<nbDimension; i++){
	    		if (maxVar < intraVariance[i]){
	    			maxVar = intraVariance[i];
	    			iVMax=i;
	    		}
	    	}
	    	
	    	c1 = new ArrayList<double[]>();
	    	c2 = new ArrayList<double[]>();
	    	
	    	avg = center[iVMax];
	    	    	
	    	for(int i=0; i<trainingData.size(); i++){
	    		double[] point  = trainingData.get(i).clone();
	    		if (point[iVMax] <= avg){
	    			c1.add(point);
	    		} else {
	    			c2.add(point);
	    		}
	    	}
	    	
	    	if(c1.size()>0) clusters.add(c1);
	    	if(c2.size()>0) clusters.add(c2);
	    	
	    	bestCluster = 0;
	    	bestInter = -1.0;
	    	for(int i=0; i<clusters.size(); i++) {
	    		intra = ArraysUtilities.getSumArray( computeIntraVariance(clusters.get(i)) );
	    		if(bestInter < intra ) {
	    			bestInter = intra;
	    			bestCluster = i;
	    		}
	    	}
	    	
	    }
	    
	    for(int c=0; c<clusters.size(); c++){
	    	resultSet.add(ArraysUtilities.getAverage(clusters.get(c), nbDimension));
	    }
	
	    return resultSet;
	}

	/**
	 * Variance Partitioning Initialization,  
	 *
	 * @param <T> type of the points to cluster
	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers
	 */
	public List<double[]> AnovaPartitioningInitialization(List<double[]> data, final int k) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	    int[] partitions = new int[numPoints];
	    Arrays.fill(partitions, 0);
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
		List<double[]> c1 = new ArrayList<double[]>(),
				c2 = new ArrayList<double[]>();
	    List < List<double[]> > clusters = new ArrayList< List<double[]> >();
	
	    // Choose one center uniformly at random from among the data points.
	    // if center is outlier, choose another
	    int nbDimension = data.get(0).length;
	    List<double[]> trainingData = new ArrayList<double[]>();
	    clusters.add(data);
	    /* compute global center */
		double[] globalCenter = ArraysUtilities.getAverage(data, nbDimension);
	    
	    int bestCluster = 0, iVMax;
		double bestInter = -1.0;
		double intra, inter, maxVar, avg, ratio;
		
		for(int c=1; c<k; c++){
	    	trainingData = clusters.get(bestCluster);
	    	
	    	clusters.remove(bestCluster);
	    	//System.out.println(trainingData.size());
	    	double[] center = ArraysUtilities.getAverage(trainingData, nbDimension);
	    	double[] intraVariance = computeIntraVariance(trainingData);
	    	double[] interVariance = computeInterVariance(globalCenter, center);
	    	maxVar = (interVariance[0]==0.0)?intraVariance[0]
	    			:intraVariance[0]/(interVariance[0]*trainingData.size());
	    	iVMax=0;
	    	for(int d=1; d<nbDimension; d++){
	    		ratio = (interVariance[d]==0.0)?intraVariance[d]
	    				:intraVariance[d]/(interVariance[d]*trainingData.size());
	    		if (maxVar < ratio){
	    			maxVar = ratio;
	    			iVMax=d;
	    		}
	    	}
	    	
	    	c1 = new ArrayList<double[]>();
	    	c2 = new ArrayList<double[]>();
	    	
	    	
	    	avg = center[iVMax];
	    		    	
	    	for(int i=0; i<trainingData.size(); i++){
	    		double[] point  = trainingData.get(i).clone();
	    		if (point[iVMax] <= avg){
	    			c1.add(point);
	    		} else {
	    			c2.add(point);
	    		}
	    	}
	    	
	    	if(c1.size()>0) clusters.add(c1);
	    	if(c2.size()>0) clusters.add(c2);
	    	
	    	bestCluster = 0;
	    	bestInter = -1.0;
	    	for(int i=0; i<clusters.size(); i++) {
	    		double[] centroid = ArraysUtilities.getAverage(clusters.get(i), nbDimension);
	    		intra = ArraysUtilities.getSumArray( computeIntraVariance(clusters.get(i)) );
	    		inter = clusters.get(i).size()*cp.distance2(centroid, globalCenter);
	    		ratio = (inter==0.0)?intra:intra/inter;
	    		if(bestInter < ratio ) {
	    			bestInter = ratio;
	    			bestCluster = i;
	    		}
	    	}
	    }
	    
	    for(int c=0; c<clusters.size(); c++){
	    	resultSet.add(ArraysUtilities.getAverage(clusters.get(c), nbDimension));
	    }
	
	    return resultSet;
	}

	/**
	 * Compute the variance of a data set on each dimension
	 * @param normalizedData data set normalized between 0 and 1
	 * @return
	 */
	private double[] computeIntraVariance(List<double[]> normalizedData) {
		
		int nbData = normalizedData.size(),
				nbDim = normalizedData.get(0).length;
		
		/* compute global center */
		double[] globalCenter = ArraysUtilities.getAverage(normalizedData, nbDim);
		
		double[] v = new double[nbDim];
		
		Arrays.fill(v, 0.0);
		
		/* compute dimension variance */
		for(int d=0; d<nbDim; d++){
			for(int i=0; i<nbData; i++){
				v[d] += Math.pow((normalizedData.get(i)[d] - globalCenter[d]), 2.0);
			}
		}
		
		return v;
	}

	/**
	 * Compute Inter-cluster Variance of a cluster on each dimension between
	 * @param globalCenter center of the data
	 * @param center center of the cluster
	 * @return
	 */
	private double[] computeInterVariance(double[] globalCenter, double[] center) {
		double[] inter = new double[globalCenter.length];	
				
		for (int d=0; d<globalCenter.length; d++){
			inter[d] = Math.pow((center[d] - globalCenter[d]), 2.0);
		}
		return inter;
	}

	/**
	
	 * Use K-means++ to choose the initial centers.
	 *
	 * @param <T> type of the points to cluster
	
	 * @param points the points to choose the initial centers from
	 * @param k the number of centers to choose
	 * @param random random generator to use
	 * @return the initial centers
	
	 */
	public List<double[]> initializationKMeansPlusPlus(List<double[]> data, final int k, double[] weights) {
	
	    // Convert to list for indexed access. Make it unmodifiable, since removal of items
	    // would screw up the logic of this method.
	    final List<double[]> pointList = Collections.unmodifiableList(new ArrayList<double[]> (data));
	
	    // The number of points in the list.
	    final int numPoints = pointList.size();
	
	    // Set the corresponding element in this array to indicate when
	    // elements of pointList are no longer available.
	    final boolean[] taken = new boolean[numPoints];
	
	    // The resulting list of initial centers.
	    final List<double[]> resultSet = new ArrayList<double[]>();
	
	    // Choose one center uniformly at random from among the data points.
	    // if center is outlier, choose another
	    int nbDimension = data.get(0).length;
	    double[] avg = ArraysUtilities.getAverage(data, nbDimension);
	    double[] sd = ArraysUtilities.getSd(data, nbDimension);
	    int tmp = 0;
	    boolean bad = true;
	    while(bad){
	    	bad=false;
	    	tmp = r.nextInt(numPoints);
	    	double[] tmpPoint = pointList.get(tmp);
	    	for(int j=0; j<nbDimension; j++){
	    		if( tmpPoint[j]>=(avg[j]+2*sd[j]) || tmpPoint[j]<=(avg[j]-2*sd[j])){
	    			bad = true;
	    		}
	    	}
	    }
	    
	    
	    final int firstPointIndex = tmp; 
	    final double[] firstPoint = pointList.get(firstPointIndex);
	
	    resultSet.add(firstPoint.clone());
	
	    // Must mark it as taken
	    taken[firstPointIndex] = true;
	
	    // To keep track of the minimum distance squared of elements of
	    // pointList to elements of resultSet.
	    final double[] minDistSquared = new double[numPoints];
	
	    // Initialize the elements.  Since the only point in resultSet is firstPoint,
	    // this is very easy.
	    for (int i = 0; i < numPoints; i++) {
	        if (i != firstPointIndex) { // That point isn't considered
	            double d = cp.weightedDistance2(firstPoint, pointList.get(i), weights);
	            minDistSquared[i] = d*d;
	        }
	    }
	
	    while (resultSet.size() < k) {
	
	        // Sum up the squared distances for the points in pointList not
	        // already taken.
	        double distSqSum = 0.0;
	
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                distSqSum += minDistSquared[i];
	            }
	        }
	
	        // Add one new data point as a center. Each point x is chosen with
	        // probability proportional to D(x)2
	        final double rn = r.nextDouble() * distSqSum;
	
	        // The index of the next point to be added to the resultSet.
	        int nextPointIndex = -1;
	
	        // Sum through the squared min distances again, stopping when
	        // sum >= r.
	        double sum = 0.0;
	        for (int i = 0; i < numPoints; i++) {
	            if (!taken[i]) {
	                sum += minDistSquared[i];
	                if (sum >= rn) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // If it's not set to >= 0, the point wasn't found in the previous
	        // for loop, probably because distances are extremely small.  Just pick
	        // the last available point.
	        if (nextPointIndex == -1) {
	            for (int i = numPoints - 1; i >= 0; i--) {
	                if (!taken[i]) {
	                    nextPointIndex = i;
	                    break;
	                }
	            }
	        }
	
	        // We found one.
	        if (nextPointIndex >= 0) {
	
	            double[] p = pointList.get(nextPointIndex);
	
	            resultSet.add(p.clone());
	
	            // Mark it as taken.
	            taken[nextPointIndex] = true;
	
	            if (resultSet.size() < k) {
	                // Now update elements of minDistSquared.  We only have to compute
	                // the distance to the new center to do this.
	                for (int j = 0; j < numPoints; j++) {
	                    // Only have to worry about the points still not taken.
	                    if (!taken[j]) {
	                        double d = cp.weightedDistance2(p, pointList.get(j), weights);
	                        double d2 = d * d;
	                        if (d2 < minDistSquared[j]) {
	                            minDistSquared[j] = d2;
	                        }
	                    }
	                }
	            }
	
	        } else {
	            // None found --
	            // Break from the while loop to prevent
	            // an infinite loop.
	            break;
	        }
	    }
	
	    return resultSet;
	}
	
	
}
