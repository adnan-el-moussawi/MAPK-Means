package Main;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;

import clustering.KMeans;
import clustering.KMeansCentroidInitialization;
import clustering.MAPKMeans;
import clustering.WKMeans;
import evaluation.ExternalEvaluation;
import readers.DoubleCSVReader;
import readers.CSVReader;
import utilities.ArraysUtilities;


/**
 * 
 * This class allows to reproduce the experiments on UCI Data presented in the section 5 of (A. EL MOUSSAWI et al., 2019).
 * <br>
 * This experiment allows to 
 * 
<ul>
	<li>
	Address the problem of tuning of the parameter &alpha;, more precisely 
	we answer the following question: is 0.5 a good default value for &alpha;?
	</li>

	<li>
	Study the impact of using preferences on the attribute preferences on 
	the clustering results. In particular, we study the sensibility of MAPK-Means 
	to the user preferences in order to answer three 	main questions: 
	(i) is there a positive impact for a relevant user-specified preferences? 
	(ii) can MAPK-means correct a bad choice of preferences?
	(iii) who better guides the exploration of the solution: user preferences or data?
	</li>

</ul>
<h3>Experimental Protocol</h3>
For each data set :
<ul>
	<li>get 100 different centers initialization using K-Means++</li>
	<li>Using two different initializations of weights (Optimal ANOVA, Bad ANOVA)</li>
	<li>for &alpha; between 0 and 1 (step = 0.1)
		<ul>
			<li>for &kappa;  between 0 and 1 (step = 0.1)</li>
				<ul>
					<li>We run MAPK-Means using the 100 initialization and for each run we compute the NMI, the objective function</li>
					<li>Over the 100 runs we compute the standard deviation and the average of NMI</li>
				</ul>
			<li>We keep the best NMI average obtained for each &alpha; value and over all &kappa; values and the corresponding 
			standard deviation and &kappa;_best</li>
		</ul>
	</li>
	<li>We keep the best NMI average obtained over all &alpha; and &kappa; values, the corresponding 
			standard deviation &alpha;_best and &kappa;_best</li>
</ul>
<br>
Finally, from those results we store in different files the following :

<ul>
	<li>A table that compares the best average of NMI values between the clusterings obtained for &alpha;=0.5 and &kappa; varying in [0; 1] vs 
	the best NMI values over all possible &alpha; and &kappa; combinations.
	(for all datasets and different initializations of weights)
	The aim is to demonstrate that 0.5 can be used as a default value of &alpha;.
	</li>
	
	<li>A table that compares the average of NMI scores while &kappa; in {0, 0.25, 0.5, 0.75, 1} 
	and using two different initializations of weights (&alpha; is set to 0.5).</li>
	
	<li>A table containing the best NMI average obtained for each &alpha; and over all &kappa; values, with 
	the corresponding standard deviation and &kappa;_best</li>
	
	<li>A table containing the results obtained for &alpha;=0.5 and &kappa; ranging in [0, 1] (step=0.05).</li>
	
</ul>

 * 
 * @author Adnan EL MOUSSAWI
 *
 */
//ParametersTuning_SensibilityToPreferencesExperiment
public class ParametersTuning_SensibilityToPreferencesExperiment {

	private static Date date;
	private static SimpleDateFormat dateFormat = new SimpleDateFormat("YYYY-MM-dd");
	private static String resultFolderPath = "../results/MAPKMeans_SensibilityToPreferencesExperiment/";
	private static String currentDate= "";
	static BufferedWriter resultFileAlphaTuningExperiment, 
	resultFilePreferencesTuningExperiment;
	
	
	public static void main(String[] args) {
		CSVReader<double[]> reader = new DoubleCSVReader();
		reader.setSep(";");
		try {
			date = new Date();
			currentDate = dateFormat.format(date);resultFolderPath+="/test "+currentDate;
			//resultFolderPath = "/home/adnan/Ictai Journal Result/alpha parameter tuning/"+dateCourant;
			
			File dir = new File(resultFolderPath);
			if (!dir.exists())
				new File(resultFolderPath).mkdirs();
			
			
			System.out.println("Debut : "+new Date());
			
			String inputPath = "./data", extension = "csv";
			String dataSets[] = new String[] {
					"Abalone",
					"Glass",
					"Iris",
					"Ionosphere",
					"Optdigits",
					"Pendigits",
					"Pgblocks",
					"Pima",
					"Vowel",
					"Wdbc",
					"Wine"
			};
			
			for(String dataset : dataSets) {
				experimentStart(reader, inputPath+"/"+dataset+"."+extension);
			}
			
			System.out.println("Fin : "+new Date()+"\n");
			

			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Initializing the protocol and start the experiments
	 * @param reader A {@link CSVReader} to read the data file
	 * @param filename the file name
	 * @throws IOException
	 */
	private static void experimentStart(CSVReader<double[]> reader, String filename) throws IOException {
		reader.open(filename);		
		List<double[]> ds = reader.read(filename);
		List<String> groundTruthPartition = reader.getLabels();
		String ShortFileName = reader.getShortName(filename);
		reader.close();
		resultFileAlphaTuningExperiment = new BufferedWriter(
				new FileWriter(resultFolderPath+"/alphaTuning - "+ShortFileName+".csv"));		
		resultFileAlphaTuningExperiment.write(";;Quality of clusterings for alpha=0.5 and all kappa;;;"
				+ "Quality of clusterings for all alpha and kappa;;;;\n");
		resultFileAlphaTuningExperiment.write("Data set;W* initialization;Avg;SD;kappa;"
				+ "Avg;SD;alpha;kappa;\n");
		
		resultFilePreferencesTuningExperiment = new BufferedWriter(
				new FileWriter(resultFolderPath+"/preferencesTuning - "+ShortFileName+".csv"));
		resultFilePreferencesTuningExperiment.write(
				";;kappa=1;;kappa=0.75;;kappa=0.5;;kappa=0.25;;kappa=0;\n");
		resultFilePreferencesTuningExperiment.write(
				"Data set;W* initialization;AVG;SD;AVG;SD;AVG;SD;AVG;SD;AVG;SD\n");
		
		System.out.println("\n"+ShortFileName);
		
		int nbDim = ds.get(0).length;
		//ds.size() / 5;
		int k = (new HashSet<>(groundTruthPartition)).size();
		
		List<double[]> normalizedData = ArraysUtilities.minmaxNormalisation(ds, nbDim); //pr.minmaxNormalisation();
		
		double[] optimalWeights = ComputeWeights.computeOptimalWeightAnova(normalizedData, groundTruthPartition, k);
		double[] badWeights = ComputeWeights.computeBadWeightAnova(normalizedData, groundTruthPartition, k);
		
		KMeansCentroidInitialization init = new KMeansCentroidInitialization();
		List<List<double[]>> centersPackOptimalAnova = new ArrayList<List<double[]>>(),
				centersPackBadAnova = new ArrayList<List<double[]>>();
		
		int nbRun = 100;
		/* get 100 different centroids initialization with weighted K-Means++ (using optimal weights), 
		 * then execute MAPK-Means */
		for(int i=0; i<nbRun; i++)
			centersPackOptimalAnova.add( init.initializationKMeansPlusPlus(normalizedData, k, optimalWeights) );
		
		resultFileAlphaTuningExperiment.write(ShortFileName+";"+"Optimal;");
		startMAPKMeans(normalizedData, groundTruthPartition, k, nbDim, 
				optimalWeights, ShortFileName, "Optimal ANOVA weights", centersPackOptimalAnova);
		System.out.println();
		
		/* get 100 different centroids initialization with weighted K-Means++ (using bad weights), 
		 * then execute MAPK-Means*/
		for(int i=0; i<nbRun; i++)
			centersPackBadAnova.add( init.initializationKMeansPlusPlus(normalizedData, k, badWeights) );
		resultFileAlphaTuningExperiment.write(ShortFileName+";"+"Bad;");
		startMAPKMeans(normalizedData, groundTruthPartition, k, nbDim, 
				badWeights, ShortFileName, "Bad ANOVA weights", centersPackBadAnova);
		System.out.println();
		
		resultFileAlphaTuningExperiment.flush();
		resultFileAlphaTuningExperiment.close();
		resultFilePreferencesTuningExperiment.close();
		
	}


	/**
	 * Execute MAPK-Means algo for all value of alpha and kappa over 100 different center initialization
	 * @param normalizedData the normalized data set
	 * @param groundTruth the ground truth partition
	 * @param k the number of clusters
	 * @param nbDim the number of dimensions
	 * @param weights
	 * @param dataSetName the name of the data set
	 * @param weightsInitializationMethod
	 * @param centersPack the List of  different centers initializations
	 * @throws IOException
	 */
	protected static void startMAPKMeans(List<double[]> normalizedData, List<String> groundTruth,
			int k, int nbDim, double[] weights, String dataSetName, 
			String weightsInitializationMethod, List<List<double[]>> centersPack) throws IOException {
		
		int nbRun = centersPack.size();
		
		//Evaluation eval = new Evaluation<>();
		List<int[]> partitions = new ArrayList<int[]>();
		int [] bestPartition = new int[normalizedData.size()];
		
		String strWeightsVectorHeader = "";
		for(int i=1; i<=weights.length; i++)
			strWeightsVectorHeader += "w_"+i+";";
		
		String strBestResultsForEachAlpha = "alpha;Best kappa;Average NMI; SD NMI",
				strDefaultAlphaAllKappaResults = "alpha;kappa;SD NMI;Average NMI;Best NMI (% Fobjectif);"
						+strWeightsVectorHeader+"\n";
		
		
		
		double[][] avgNMI = new double[101][101],
				sdNMI = new double[101][101];
				
		double maxNMIAvgCurrentAlpha = 0.0,
				bestKappaValueForCurrentAlpha = 0.0, alpha = 0.0, kappa, minFobjectiveForEachKappa,
				maxNMIAvgOverAllAlpha = 0.0, bestAlpha = 0.0, bestKappaOverAllAlpha = 0.0;
	
		double[] bestLearnedWeightsForEachKappa = new double[weights.length];
		int indexBestKappa;
				
		MAPKMeans algo;
		
		
		/* the confidence \kappa \in [0, 1[, with step=0.05 */
		for(int i = 0; i < 100; i=i+5) {
			
			 maxNMIAvgCurrentAlpha = -2.0;
			 alpha = ((double)i)/100;
			if( i%5 == 0 )
				System.out.print(i+"% ");
			/* the confidence \kappa \in [0, 1], with step=0.05 */
			for(int j=0; j <= 100; j=j+5){
				kappa = ((double)j)/100;				
				
				double[] nmi = new double[nbRun];
				double[] fobjective = new double[nbRun];
				int [] tmp_partition;
				minFobjectiveForEachKappa = Double.MAX_VALUE;
				
				/* execute MAPK-Means nbRun executions using different initial centers 
				 * from pre-computed list (centersPacks) */
				for(int run = 0; run < nbRun; run++){
					algo = new MAPKMeans();
					algo.cluster(normalizedData, centersPack.get(run), k, 100, weights, alpha, kappa);
					fobjective[run] = algo.getIntraClusterDistance();
					
					tmp_partition = algo.getPartition();
					nmi[run] = ExternalEvaluation.entropyNMISun2010(tmp_partition, groundTruth);
					
					/* get the learned weights vector of the best solution minimizing the objective function */
					if(fobjective[run] < minFobjectiveForEachKappa){
						bestLearnedWeightsForEachKappa = algo.getLearningVector();					
					}
				}
				
				/* compute the average and the standard deviation of the NMI scores over all runs */
				avgNMI[i][j] = ArraysUtilities.getAverageArray(nmi);
				sdNMI[i][j] = ArraysUtilities.getSdArray(nmi);
				
				partitions.add(bestPartition.clone());
				
				/* save the computed NMI values for \alpha=0.5 and each \kappa */
				if(i == 50){
					strDefaultAlphaAllKappaResults += alpha+";"
							+kappa+";"
							+sdNMI[i][j]+";"
							+avgNMI[i][j]+";"
							+avgNMI[i][j]+";"
							+toStringArray(bestLearnedWeightsForEachKappa)
							+"\n";
				}

				
				
				/* choose the best NMI average over all \kappa values and the corresponding best \kappa, 
				 * choose the largest \kappa if many best values exist */
				if(maxNMIAvgCurrentAlpha <= avgNMI[i][j]){
					maxNMIAvgCurrentAlpha = avgNMI[i][j];
					bestKappaValueForCurrentAlpha = kappa;
					
				}
				
				// To do the same but while choosing the best smallest kappa
				/*if(maxAVGCurrentAlpha < avgNMI[i][j]){
					maxAVGCurrentAlpha = avgNMI[i][j];
					bestKappaValueForCurrentAlpha = kappa;
					
				}*/
				
			}//end for all kappa
			
			
			indexBestKappa = (int) (bestKappaValueForCurrentAlpha*100);
			
			/* save the best results for each \alpha */
			strBestResultsForEachAlpha += 
					alpha+";"
					+ bestKappaValueForCurrentAlpha+";"
					+ avgNMI[i][indexBestKappa]+";"
					+ sdNMI[i][indexBestKappa]+";"
					+ "\n";
			/*
			 * get the best clustering for all alpha and kappa (i.e. maximizing the average of NMI),
			 * 
			 */
			if(maxNMIAvgOverAllAlpha <= maxNMIAvgCurrentAlpha){//choose the smallest best alpha
				//save the best NMI average
				maxNMIAvgOverAllAlpha = maxNMIAvgCurrentAlpha;
				// save the corresponding best alpha value
				bestAlpha = alpha;
				// save the corresponding best kappa value
				bestKappaOverAllAlpha = bestKappaValueForCurrentAlpha; // save its c
				
			}
			
			/*
			 * get the best clustering for alpha=0.5 and all kappa
			 */
			if(i == 50){//choose the smallest best alpha
				resultFileAlphaTuningExperiment.write(
						avgNMI[i][indexBestKappa]+";"
						+ sdNMI[i][indexBestKappa]+";"
						+ bestKappaValueForCurrentAlpha +";");
			}
			//System.out.println();			
		}//end for all alpha
		int indexBestAlpha = (int) (bestAlpha*100);
		int indexKappaAllAlpha = (int) (bestKappaOverAllAlpha*100);
		
		/* save the results for alpha tuning of this data set to the global results file (for all datasets) */
		resultFileAlphaTuningExperiment.write(
				avgNMI[indexBestAlpha][indexKappaAllAlpha]+";"
				+ sdNMI[indexBestAlpha][indexKappaAllAlpha]+";"
				+ bestAlpha+";"
				+ bestKappaOverAllAlpha+";");
		resultFileAlphaTuningExperiment.flush();
		
		/* save the results for preferences tuning of this data set to the global results file (for all datasets) */
		resultFilePreferencesTuningExperiment.write(dataSetName+";"
				+weightsInitializationMethod+";"
				+ avgNMI[50][100]+";"
				+ sdNMI[50][100]+";"
				+ avgNMI[50][75]+";"
				+ sdNMI[50][75]+";"
				+ avgNMI[50][50]+";"
				+ sdNMI[50][50]+";"
				+ avgNMI[50][25]+";"
				+ sdNMI[50][25]+";"
				+ avgNMI[50][0]+";"
				+ sdNMI [50][0]+"\n");
		resultFilePreferencesTuningExperiment.flush();
		
		/* save some particular results for the current dataset and the current weights iniialization method */
		BufferedWriter file1BestResultsForAllAlphaValues, 
		file2ResultsForDefaultAlphaAndAllKappaValues;
		
		String repertoire = resultFolderPath+"/"+dataSetName+"/MAPK-Means "+weightsInitializationMethod;
		File dir = new File(repertoire);
		if (!dir.exists())
			new File(repertoire).mkdirs();
		
		// store the best results for all alpha values
		file1BestResultsForAllAlphaValues = new BufferedWriter(
				new FileWriter(repertoire+"/Best results for All alpha values - "+currentDate+".csv"));
		file1BestResultsForAllAlphaValues.write(strBestResultsForEachAlpha+"\n");
		file1BestResultsForAllAlphaValues.flush();
		file1BestResultsForAllAlphaValues.close();
		
		// store the best results for alpha=0.5 and all kappa values
		file2ResultsForDefaultAlphaAndAllKappaValues = new BufferedWriter(
				new FileWriter(repertoire+"/Results for alpha set to 0.5 - "+currentDate+".csv"));
		file2ResultsForDefaultAlphaAndAllKappaValues.write(strDefaultAlphaAllKappaResults+"\n");
		file2ResultsForDefaultAlphaAndAllKappaValues.flush();
		file2ResultsForDefaultAlphaAndAllKappaValues.close();
	}



	/* *************************************** */
	static String toStringArray(int[] a){
		String str = "";
		for(int i=0; i<a.length; i++) str += a[i]+";";
		return str;
	}
	
	static String toStringArray(double[] a){
		String str = "";
		for(int i=0; i<a.length; i++) str += a[i]+";";
		return str;
	}
	
	static String toStringArray(float[] a){
		String str = "";
		for(int i=0; i<a.length; i++) str += a[i]+";";
		return str;
	}
	
	static String toStringArray2(double[] a){
		String str = "[";
		for(int i=0; i<a.length; i++) str += a[i]+",";
		str = str.substring(0, str.length()-1);
		return str+"]";
	}
}
