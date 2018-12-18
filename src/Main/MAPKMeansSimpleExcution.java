package Main;
import java.util.Date;
import java.util.List;

import clustering.MAPKMeans;
import readers.DoubleCSVReader;
import readers.CSVReader;
import utilities.ArraysUtilities;

public class MAPKMeansSimpleExcution {

	private static CSVReader<double[]>  reader;

	public static void main(String[] args){		
		reader = new DoubleCSVReader();
		reader.setSep(";");
		
		
		
		System.out.println("Initialisation K-Means++...");
		System.out.println("Debut : "+new Date());
		
		// reading the data
		String filename = "./data/ionosphere_o.csv";
		reader.open(filename );		
		List<double[]> ds = reader.read(filename);
		reader.close();
		
		//
		int nbDim = ds.get(0).length;
		int k = 3;
		double kappa=1;
		
		// normalize the data between 0 and 1
		List<double[]> normalizedData = ArraysUtilities.minmaxNormalisation(ds, nbDim); //pr.minmaxNormalisation();
		
		// running MAPK-Means with default centers initialization (KMeans++)
		MAPKMeans algo = new MAPKMeans();
		//algo.cluster(normalizedData, k, 100, poids, 0.5, kappa); // for this method poids should be set 
		algo.clusterUniformInit(normalizedData, k, 100, null, 0.5, kappa);
		
		
		// running MAPK-Means with personalized centers initialization (KMeans++)
		//List<double[]> centres = initialization.initializationKMeansPlusPlus(data, kint, poids); // KMeans++ initialization
		//List<double[]> centres = initialization.AnovaPartitioningInitialization(data, kint); // ANOVA initialization
		//algo.cluster(normalizedData, centers, k, 100, poids, 0.5, kappa);
		//algo.clusterUniformInit(normalizedData, centers, k, 100, null, 0.5, kappa);
		
		
		
		// getting the results
		int[] partition = algo.getPartition(); //data partition (vector of clusters labels)
		double dIntraCluster = algo.getIntraClusterDistance(); // total intra-clusters distances
		double fobjective = algo.getObjectiveFunction(); // objective function value
		double[] learnedWeigths = algo.getLearningVector(); // final vector of weights (computed by the algorithm)
		
		System.out.println("final weights : "+toStringArray(learnedWeigths));
		
		System.out.println("Fin : "+new Date()+"\n");
		
	}
	
	static String toStringArray(double[] a){
		String str = "";
		for(int i=0; i<a.length; i++) str += a[i]+";";
		return str;
	}
}
