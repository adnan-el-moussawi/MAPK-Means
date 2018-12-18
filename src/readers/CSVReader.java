package readers;

import java.util.List;
/**
 * Read a CVS file with ONLY numeric data
 * @author Adnan
 *
 * @param <D> a primitive Java data numeric type
 */
public interface CSVReader<D> {
	/**
	 * Open a CSV file, read, parse then return the data
	 * @param filePath The file path
	 * @return The numeric data
	 */
	public List<D> read(String filePath);
	
	/**
	 * Open a CSV file, read then return the data
	 * @param filePath The file path
	 */
	public List<D> read(String filePath, boolean hasLables);
	
	/**
	 * Open a CSV file without reading the data
	 * @param filePath The file path
	 */
	public void open(String filePath);
	
	/**
	 * Read the next n points (n lines) of the data file, if exists
	 * @param n The number of lines to read
	 */
	public void next(int n);

	/**
	 * Reset the reader and the counter (to read a new file)
	 */
	public void reset();	
	
	/**
	 * close a data file
	 */
	public void close();

	/**
	 * Return the number of current redden data lines
	 * @return 
	 */
	public int count();

	/**
	 * Get the total number of data points (number lines) after using {@link IReader.open()} or {@link IReader.read()} 
	 * @return
	 */
	public int size();
	
	/**
	 * Get the data set
	 * @return A list of numeric data for example List<double[]>
	 */
	public List<D> getDataSet();

	/**
	 * Get the column representing the clusters labels
	 * @return A String list of clusters labels exiting in the last column  of the CSV file
	 */
	public List<String> getLabels();

	/**
	 * Get the short name of the data file
	 * @param filePath The file path
	 * @return The short file name
	 */
	public default String getShortName(String filePath) {
		int idx = filePath.replaceAll("\\\\", "/").lastIndexOf("/");
		int idx2 = filePath.lastIndexOf(".");
		return idx >= 0 ? filePath.substring(idx + 1, idx2) : filePath;
	}

	/**
	 * Define the CSV data separator
	 * @param sep The data separator (";" by default)
	 */
	public void setSep(String sep);;
}
