package readers;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * Read CVS file with numeric data with (or without) clusters labels
 * <br>data format : 
 * <br>first line contains the name of attributes than each line contains numeric values and class name separated by ";" 
 * <br>a1;a2;a3;class_labels
 * <br>219;22;21.5;class1
 * <br>211;25;21.6;class1
 * <br>221;21;22.8;class1
 * <br>125;32;1.15;class2
 * <br>126;34;1.85;class2
 * @author Adnan EL MOUSSAWI
 */
public class DoubleCSVReader implements CSVReader<double[]>{
	/**
	 * Data separator (";" by default) 
	 */
	private String sep = " ";
	
	/**
	 * {@link List} of double array to store the data
	 */
	List<double[]> data;
	
	/**
	 * {@link List} of {@link String} to store clusters labels if exists
	 */
	List<String> labels;

	/**
	 * An array of string to the name of attributes (the file header)
	 */
	String[] attributeNames;
	
	/**
	 * A {@link BufferedReader} to read the data file
	 */
	BufferedReader br;

	/**
	 * A counter to store the number of redden lines
	 */
	int count;
	
	/**
	 * A string to store the last redden line
	 */
	String line;
	
	/**
	 * File path
	 */
	String filePath;

	public DoubleCSVReader() {}
	
	/**
	 * @param sep data separator, example ";" or ","
	 */
	public DoubleCSVReader(String sep){
		this.setSep(sep);
	}
	
	/**
	 * Open a CSV file without reading the data
	 * @param filePath The file path
	 */
	public void open(String filePath) {
		try{
			this.filePath = filePath;
			this.count = 0;
			this.data = new ArrayList<double[]>();
			this.labels = new ArrayList<String>();
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
			attributeNames = br.readLine().split(sep);
			while ((line = br.readLine())!=null){
	            count ++;
	        }
	
			br.close();
			// Remise du flux au debut
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
	
		} catch (Exception e) {e.printStackTrace();}
	}
	
	/**
	 * Open a CSV file, read then return the data
	 * @param filePath The file path
	 */
	public List<double[]> read(String filePath) {
		try{
			this.count = 0;
			this.data = new ArrayList<double[]>();
			this.labels = new ArrayList<String>();
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
			br.readLine();
			while ((line = br.readLine())!=null){
				parseRecord(line);
	            count ++;
	        }
			return this.data;
		} catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Parse a CSV record and save the data
	 * @param line A CSV record (line)
	 */
	private void parseRecord(String line){
		String[] tmp = line.split(sep);
		int i;
		double[] vect = new double[tmp.length - 1];
		for (i = 0; i < tmp.length - 1; i++){
			vect[i] = Double.parseDouble(tmp[i]);
		}
		
		this.data.add((double[]) vect);
		this.labels.add(tmp[tmp.length - 1]);
	}
	
	/**
	 * Open a CSV file, read then return the data
	 * @param filePath The file path
	 */
	public List<double[]> read(String filePath, boolean hasLables) {
		if(hasLables) {
			return read(filePath);
		}
		else {
			try{
				this.count = 0;
				this.data = new ArrayList<double[]>();
				this.labels = null;
				br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
				br.readLine();
				while ((line = br.readLine())!=null){
					parseRecordWithoutLabel(line);
		            count ++;
		        }
				return this.data;
			} catch(Exception e){
				e.printStackTrace();
				return null;
			}
		}
	}

	/**
	 * Parse a CSV record and save the data
	 * @param line A CSV record (line)
	 */
	private void parseRecordWithoutLabel(String line){
		String[] tmp = line.split(sep);
		int i;
		double[] vect = new double[tmp.length];
		for (i = 0; i < tmp.length; i++){
			vect[i] = Double.parseDouble(tmp[i]);
		}
		
		this.data.add((double[]) vect);
	}
	
	/**
	 * Read the next n points (n lines) of the data file, if exists
	 * @param n The number of lines to read
	 */
	public void next(int n) {
		this.data.clear();
		this.labels.clear();
		try {
			int nb = 0; // Nombre de lignes lues
			boolean eof = false;
			while (nb < n && !eof){
				line = this.br.readLine();
				if (line !=null){
					parseRecord(line);
					count ++;
					nb++;
				} else {
					eof = true;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void close() {
		try{
			if (br != null) br.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**
	 * Reset the reader and the counter (to read a new file)
	 */
	public void reset() {
		try {
			if (br != null) {
				br.close();
			}
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public int size() {
		return getDataSet().size();
	}
	
	/**
	 * Return the number of current redden data lines
	 * @return 
	 */
	public int count() {
		return this.count;
	}

	/**
	 * Get the numeric data
	 */
	public List<double[]> getDataSet() {
		return this.data;
	}
	
	/**
	 * Get clusters labels
	 */
	public List<String> getLabels() {
		return this.labels;
	}
	
	/**
	 * Set the CSV data separator
	 * @param sep Data separator (";", ",", "\t" or " ")
	 */
	public void setSep(String sep){
		this.sep = sep;
	}
}
