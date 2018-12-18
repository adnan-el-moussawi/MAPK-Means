/**
 * 
 */
package utilities;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This class provides many methods for math operators on numeric arrays and list of of arrays ({@link List<double[]>}
 * @author ADNAN EL MOUSSAWI
 *
 */
public class ArraysUtilities implements Serializable {
	/**
	 * serialization UID
	 */
	private static final long serialVersionUID = -9064395469452434284L;

	
	/**
	 * Get the sum of data on all dimensions
	 * @param ds The data set, a list of data points 
	 * @param nbDimension The number of dimensions
	 * @return An array of the sum of data on each dimension
	 */
	public static double[] getSomme(List<double[]> ds, int nbDimension){
		
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			return array;
		}
		else if (ds.size()==1){
			return (ds.get(0).clone());
		}
		else{
			double[] sum = new double[nbDimension];
			Arrays.fill(sum, 0.0);
			double[] array;
			for(int i=0; i<ds.size(); i++){
				array = ds.get(i);
				for(int j=0; j<nbDimension; j++){
					sum[j]+=array[j];
				}
			}
			return sum;
		}
	}

	/**
	 * Get the minimum values for all dimensions for a data set
	 * @param ds The data set, a list of data points 
	 * @param nbDimension The number of dimensions
	 * @return An array of the minimum values for each dimension
	 */
	public static double[] getMin(List<double[]> ds, int nbDimension){
				
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			return array;
		}
		else{

			double[] min = ds.get(0).clone();
			double[] array = new double[nbDimension];
			for(int i=1; i<ds.size(); i++){
				array = ds.get(i);
				for(int j=0; j<nbDimension; j++){
					if(min[j] > array[j]){
						min[j]=array[j];
					}
				}
			}//System.out.println(Arrays.toString(min));
			return (min.clone());
		}		
	}
	
	/**
	 * Get the maximum values for all dimensions for a data set
	 * @param ds The data set, a list of data points 
	 * @param nbDimension The number of dimensions
	 * @return An array of the maximum values for each dimension
	 */
	public static double[] getMax(List<double[]> ds, int nbDimension){
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			return array;
		}
		else{
			double[] max = ds.get(0).clone();
			double[] array = new double[nbDimension];
			for(int i=1; i<ds.size(); i++){
				array = ds.get(i);
				for(int j=0; j<nbDimension; j++){
					if(max[j]<array[j])
						max[j]=array[j];
				}
			}
			//System.out.println(Arrays.toString(max));
			return max.clone();
		}
	}

	
	/**
	 * Get the average values for all dimensions for a data set
	 * @param ds The data set, a list of data points 
	 * @param nbDimension The number of dimensions
	 * @return An array of the average values of each dimension
	 */
	public static double[] getAverage(List<double[]> ds, int nbDimension){
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			return array;
		}
		else if (ds.size()==1){
			return (ds.get(0).clone());
		}
		else{
			double[] sum = new double[nbDimension];
			Arrays.fill(sum, 0.0);
			double[] array;
			for(int i=0; i<ds.size(); i++){
				array = ds.get(i);
				for(int j=0; j<nbDimension; j++){
					sum[j]+=array[j];
				}
			}
			for(int j=0; j<nbDimension; j++){
				sum[j] = sum[j]/ds.size();
			}
			return sum;
		}
	}

	
	/**
	 * Get the average value of a particular dimension for a data set
	 * @param ds The data set, a list of data points
	 * @param index The index of the dimension
	 * @return An array of the average value of the dimension 
	 */
	public static double getDimensionAverage(List<double[]> ds, int index){
		if(ds.size()==0){
			
			return 0.0;
		}
		else if (ds.size()==1){
			return (ds.get(0)[index]);
		}
		else{
			double sum = 0.0;
			for(int i=0; i<ds.size(); i++){
				sum += ds.get(i)[index];
			}
			sum = sum/ds.size();
			return sum;
		}
	}
	
	/**
	 * Get the standard deviation values for all dimensions of a data set
	 * @param ds The data set, a list of data points
	 * @param nbDimension The number of dimensions
	 * @return An array of the standard deviation values for each dimension
	 */
	public static double[] getSd(List<double[]> ds, int nbDimension){
		double[] average = getAverage(ds, nbDimension);
		return getSd(ds, average, nbDimension);
	}

	/**
	 * Get the standard deviation values for all dimensions of a data set
	 * @param ds The data set, a list of data points
	 * @param average array of the average values for each dimension of The data set
	 * @param nbDimension The number of dimensions
	 * @return An array of the standard deviation values for each dimension
	 */
	public static double[] getSd(List<double[]> ds, double[] average, int nbDimension){
		double[] sd = new double[nbDimension];
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			return array;
		}
		else if (ds.size()==1){
			return (ds.get(0).clone());
		}		
		else{
			Arrays.fill(sd, 0.0);
			double[] array;
			for(int i=0; i<ds.size(); i++){
				array = ds.get(i);
				for(int j=0; j<nbDimension; j++){
					sd[j] += Math.pow((array[j]-average[j]), 2);
				}
			}
			for(int j=0; j<nbDimension; j++){
				sd[j] = (double) Math.sqrt(sd[j]/(ds.size()));
			}
			//System.out.println("sd: "+ Arrays.toString(sd));
			return sd;
		}
	}
	

	
	/**
	 * Get the standard deviation value for a dimension of The data set
	 * @param ds The data set, a list of data points
	 * @param nbDimension The number of dimensions
	 * @param index The index of the dimension
	 * @return An array of the standard deviation values for each dimension
	 */
	public static double getSd(List<double[]> ds, int nbDimension, int index){
		double average = getDimensionAverage(ds, index);
		return getSd(ds, average, nbDimension, index);
	}

	/**
	 * Get the standard deviation value for a dimensions of The data set
	 * @param ds The data set, a list of data points
	 * @param average An array of the average values for each dimension of The data set
	 * @param nbDimension The number of dimensions
	 * @param index The index of the dimension
	 * @return An array of the standard deviation values each all dimension
	 */
	public static double getSd(List<double[]> ds, double average, int nbDimension, int index){
		double sd = 0.0;
		if(ds.size()==0){
			return 0.0;
		}
		else if (ds.size()==1){
			return (ds.get(0)[index]);
		}		
		else{
			
			for(int i=0; i<ds.size(); i++){
				sd += Math.pow((ds.get(i)[index]-average), 2);
			}
			
			sd = (double) Math.sqrt(sd/(ds.size()));
			//System.out.println("sd: "+ Arrays.toString(sd));
			return sd;
		}
	}

	
	/**
	 * Normalize a data set between 0 and 1 using MinMax normalization 
	 * @param ds The data set, a list of data points
	 * @param nbDimension The number of dimensions
	 * @return Normalized data set between 0 and 1
	 */
	public static List<double[]> minmaxNormalisation(List<double[]> ds, int nbDimension){
		List<double[]> nDs = new ArrayList<double[]>();
		if(ds.size()>1){		
			double[] min = getMin(ds, nbDimension);
			double[] max = getMax(ds, nbDimension);
			//System.out.println(Arrays.toString(min));
			//System.out.println(Arrays.toString(max));
			
			double[] diff = new double[nbDimension];
			for(int j=0; j<nbDimension; j++){
				diff[j] = max[j]-min[j];
			}
			
			double[] array = new double[nbDimension];
			for(int i=0; i<ds.size(); i++){
				for(int j=0; j<nbDimension; j++){
					array[j] = (ds.get(i)[j]-min[j])/diff[j];
				}
				//System.out.println(Arrays.toString(ds.get(i)));
				nDs.add(array.clone());
			}
		}
		else if (ds.size()==1){
			nDs.add(ds.get(0).clone());
		}
		else{
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			nDs.add(array.clone());
		}
		return nDs;
	}
	
	/**
	 * Normalize a data set between 0 and 1 using MinMax normalization and pre-computed MinMax arrays 
	 * @param ds The data set, a list of data points
	 * @param min An array of the minimum values for all dimensions
	 * @param max An array of the maximum values for all dimensions
	 * @param nbDimension The number of dimensions
	 * @return Normalized data set between 0 and 1
	 */
	public static List<double[]> minmaxNormalisation(List<double[]> ds, double[] min, double[] max, int nbDimension){
		List<double[]> nDs = new ArrayList<double[]>();
		if(ds.size()>1){
			double[] diff = new double[nbDimension];
			for(int j=0; j<nbDimension; j++){
				diff[j] = max[j]-min[j];
			}
			
			double[] array = new double[nbDimension];
			for(int i=0; i<ds.size(); i++){
				for(int j=0; j<nbDimension; j++){
					array[j] = (ds.get(i)[j]-min[j])/diff[j];
				}
				//System.out.println(Arrays.toString(ds.get(i)));
				nDs.add(array.clone());
			}
		}
		else if (ds.size()==1){
			nDs.add(ds.get(0).clone());
		}
		else{
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			nDs.add(array.clone());
		}
		return nDs;
	}
	
	/**
	 * Normalize a data set between using the ZScore normalization method
	 * @param ds The data set, a list of data points
	 * @param nbDimension The number of dimensions
	 * @return Normalized data set
	 */
	public static List<double[]> zscoreNormalisation(List<double[]> ds, int nbDimension){
		List<double[]> nDs = new ArrayList<double[]>();
		if(ds.size()==0){
			double[] array = new double[nbDimension];
			Arrays.fill(array, 0.0);
			nDs.add(array.clone());
		}
		else if (ds.size()==1){
			nDs.add(ds.get(0).clone());
		}
		else{
			double[] average = getAverage(ds, nbDimension);
			double[] sd = getSd(ds, average, nbDimension);
			double[] array;
			for(int i=0; i<ds.size(); i++){
				array = new double[nbDimension];
				for(int j=0; j<nbDimension; j++){
					array[j] = (ds.get(i)[j]-average[j])/sd[j];
				}
				//System.out.println(Arrays.toString(ds.get(i)));
				nDs.add(array);
			}
		}
		return nDs;
	}
	
	
	/**
	 * Get the maximum value of an {@link int} array
	 * @param a An array of the int values
	 * @return Maximum value of the array
	 */
	public static int getMaxArray(int[] a){
		int mx = a[0];
		for (int i = 1; i < a.length; i++){
			if (a[i] > mx) mx = a[i];
		}
		return mx;
	}

	/**
	 * Get the minimum value of an {@link int} array
	 * @param a An array of the int values
	 * @return Minimum value of the array
	 */
	public static int getMinArray(int[] a){
		int min = a[0];
		for (int i = 1; i < a.length; i++){
			if (a[i] < min) min = a[i];
		}
		return min;
	}
	
	/**
	 * Get the maximum value of a {@link double} array
	 * @param a An array of the double values
	 * @return Maximum value of the array
	 */
	public static double getMaxArray(double[] a) {
		double max = a[0];
		for(int d=1; d<a.length; d++){
			if(max<a[d]){
				max=a[d];
			}
		}
		return max;
	}
	
	/**
	 * Get the minimum value of a {@link double} array
	 * @param a An array of double values
	 * @return Minimum value of the array
	 */
	public static double getMinArray(double[] a) {
		double min = a[0];
		for(int d=1; d<a.length; d++){
			if(min>a[d]){
				min=a[d];
			}
		}
		return min;
	}
	
	/**
	 * Get the sum of the elements of an array
	 * @param a an array of int values
	 * @return The sum of the array element values
	 */
	public static int getSumArray(int[] a){
		int sum = a[0];
		for (int i = 1; i < a.length; i++){
			sum += a[i];
		}
		return sum;
	}
	
	/**
	 * Get the sum of the elements of an array
	 * @param a An array of double values
	 * @return The sum of the array element values
	 */
	public static double getSumArray(double[] a){
		double sum=a[0];
		for(int d=1; d<a.length; d++){
			sum += a[d];
		}
		return sum;
	}
	
	/**
	 * Get the average of an array
	 * @param a An array of int values
	 * @return The average of the array elements
	 */
	public static double getAverageArray(int[] a){
		return (getSumArray(a)/a.length);
	}
	
	/**
	 * Get the average of an array
	 * @param a An array of the double values
	 * @return The average of the array elements
	 */
	public static double getAverageArray(double[] a){
		return (getSumArray(a)/a.length);
	}
	
	/**
	 * Get the standard deviation of an array
	 * @param a an array of int values
	 * @return The standard deviation of the array elements
	 */
	public static double getSdArray(int[] a){
		double average = getAverageArray(a);
		if(average==0.0) return 0.0;
		else{
			double sd =0.0;
			for(int j=0; j<a.length; j++){
				sd += Math.pow((a[j]-average), 2);
			}
			sd = (double) Math.sqrt(sd/(a.length));			
			return sd;
		}
	}
	
	/**
	 * Get the standard deviation of an array
	 * @param a An array of double values
	 * @return The standard deviation of the array elements
	 */
	public static double getSdArray(double[] a){
		double average = getAverageArray(a);
		if(average==0.0) return 0.0;
		else{
			double sd =0.0;
			for(int j=0; j<a.length; j++){
				sd += Math.pow((a[j]-average), 2);
			}
			
			sd = (double) Math.sqrt(sd/(a.length));
			
			return sd;
		}
	}
	
	/**
	 * Get the sum of two arrays
	 * @param a1 An array of double values
	 * @param a2 An array of double values
	 * @return The sum of two values
	 */
	public static double[] sumTwoArray(double[] a1, double[] a2){
		double[] sum = a1.clone();
		for(int d=0; d<a1.length; d++){
			sum[d] += a2[d];
		}
		return sum;
	}
	
	/**
	 * Get the subtraction of two arrays
	 * @param a1 An array of double values
	 * @param a2 An array of double values
	 * @return The difference between two arrays
	 */
	public static double[] subTwoArray(double[] a1, double[] a2){
		double[] sum = a1.clone();
		for(int d=0; d<a1.length; d++){
			sum[d] -= a2[d];
		}
		return sum;
	}
	
	/**
	 * Get the product of two arrays with same length : element by element, a1[i]*a2[i].
	 * @param a1 An array of double values
	 * @param a2 An array of double values
	 * @return The product of two arrays : element by element, a1[i]*a2[i] Example
	 * <br>a1[] = new double[]{1, 2, 3}
	 * <br>a2[] = new double[]{4, 5, 6}
	 * <br>productTwoArray(a1,a2) : new double[]{4, 10, 18}
	 */
	public static double[] productTwoArray(double[] a1, double[] a2){
		double[] sum = new double[a1.length];
		for(int d=0; d<a1.length; d++){
			sum[d] = a1[d]*a2[d];
		}
		return sum;
	}
	
	/**
	 * Get the division of the first array by the second one : element by element, a1[i]/a2[i].
	 * <br> Arrays should have same length, not null
	 * <br> If an element from a2 equal to zero, the result of the corresponding division is set to 0 
	 * @param a1 An array of double values
	 * @param a2 An array of double values
	 * @return The product of two arrays : element by element, a1[i]*a2[i] Example
	 * <br>a1[] = new double[]{1, 2, 3}
	 * <br>a2[] = new double[]{4, 5, 6}
	 * <br>divideTwoArray(a1,a2) : new double[]{0.25, 0.4, 0.5}
	 */
	public static double[] divideTwoArray(double[] a1, double[] a2){
		double[] sum = new double[a1.length];
		for(int d=0; d<a1.length; d++){
			if(a2[d]==0) sum[d]=0;
			else sum[d] = a1[d]/a2[d];
		}
		return sum;
	}
	
	/**
	 * Returns a vector with each element equal to corresponding element of 
	 * the first argument raised to the power of the second argument. 
	 * (based on {@link Math pow(a, b)}
	 * <br> Arrays should have same length, not null
	 * @param a An array of double values
	 * @param b The exponent
	 * @return A vector with each element equal to a[i]^b. Example :
	 * <br>a[] = new double[]{1, 2, 3}; b=2.0
	 * <br>powerArray(a,b) : new double[]{1, 4, 9}
	 */
	public static double[] powerArray(double[] a, double b){
		double[] sum = new double[a.length];
		for(int d=0; d<a.length; d++){
			sum[d] = Math.pow(a[d], b);
		}
		return sum;
	}
	
	/**
	 * Get the difference between two arrays (in percentage)
	 * @param a1 The first array of double values
	 * @param a2 The second array of double values
	 * @return Percentage difference between two arrays (for each dimension)
	 */
	public static double[] percentDiffBetweenTwoArray(double[] a1, double[] a2){
		double[] sum = a1.clone();
		for(int d=0; d<a1.length; d++){
			sum[d] = 100*(a2[d] - a1[d])/(a1[d]);
		}
		return sum;
	}
	
	/**
	 * Normalize array elements such their sum is equal to 1
	 * @param a An array of double values
	 * @return A normalized double array
	 */
	public static double[] normalizeArray(double[] a) {
		double[] an = new double[a.length];
		double sum = getSumArray(a);
		for (int i = 0; i < a.length; i++)
			an[i] = a[i] / sum;
		return an.clone();
	}
	
	/**
	 * Normalize array elements such their sum is equal to 1
	 * @param a An array of {@link String} (get double value from String)
	 * @return A normalized double array
	 */
	public static double[] normalizeArray(String[] a) {
		double[] tmp = toDoubleArray(a);
		return normalizeArray(tmp);
	}
	
	/**
	 * Get double array values from a {@link String} array 
	 * @param a An array of {@link String}
	 * @return An array of double values
	 */
	private static double[] toDoubleArray(String[] a) {
		double[] a_double = new double[a.length];
		for (int i = 0; i < a.length; i++)
			a_double[i] = Double.parseDouble(a[i]);
		return a_double.clone();
	}

	/**
	 * Sort a {@link double} array in ascendant order
	 * @param vector The array to sort
	 * @return Return a sorted array, in ascendant order
	 */
	public static double[] sortAsc(double vector[]) {
		double[] sortedTable = vector.clone();
		int vectorLength = vector.length;
		double tempValue = 0.0;
		boolean permut;
		
		for (int j=0; j<vectorLength; j++)
		/**
		 * sort the table ascending
		 */
		do {
			permut = false;
			for (int i = 0; i < vectorLength - 1; i++) {
				// check if two adjacent values are sorted 
				if (sortedTable[i] > sortedTable[i + 1]) {
					// if not swap position
					tempValue = sortedTable[i];
					sortedTable[i] = sortedTable[i + 1];
					sortedTable[i + 1] = tempValue;
					permut = true;
				}
				
			}
		} while (permut);
		
		return sortedTable;
	}

	/**
	 * Sort a {@link double} array in descendant order
	 * @param vector The array to sort
	 * @return Return a sorted array, in descendant order
	 */
	public static double[] sortDesc(double vector[]) {
		double[] sortedTable = vector.clone();
		int vectorLength = vector.length;
		double tempValue = 0.0;
		boolean permut;
		
		for (int j=0; j<vectorLength; j++)
		/**
		 * sort the table descending
		 */
		do {
			permut = false;
			for (int i = 0; i < vectorLength - 1; i++) {
				// check if two adjacent values are sorted 
				if (sortedTable[i] < sortedTable[i + 1]) {
					// if not swap position
					tempValue = sortedTable[i];
					sortedTable[i] = sortedTable[i + 1];
					sortedTable[i + 1] = tempValue;
					permut = true;
				}
				
			}
		} while (permut);
		
		return sortedTable;
	}
}
