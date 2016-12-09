package vis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ProblemSpec {
	/** True iff a problem is currently loaded */
	private boolean problemLoaded = false;
	/** True iff a solution is currently loaded */
	private boolean solutionLoaded = false;
	public int numObjects;
	public int numTimeSteps;
	private int timeC;
	private int objC;
	public Map<Integer, List<ObjConfig>> timesteps = new HashMap<Integer, List<ObjConfig>>();



	/**
	 * Loads a problem from a problem text file.
	 * 
	 * @param filename
	 *            the path of the text file to load.
	 * @throws IOException
	 *             if the text file doesn't exist or doesn't meet the assignment
	 *             specifications.
	 */
	public void loadProblem(String filename) throws IOException {
		problemLoaded = false;
		solutionLoaded = false;
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));

			//Read in number of objects and timesteps
			String line = br.readLine();
			//numV = Integer.parseInt(line);
			String delim = " ";
			String[] tokens = line.split(delim);
			int tokenCount = tokens.length;
			if (tokenCount != 2){
				System.out.print("Invalid environment file\n");
				System.exit(1);
			}
			numObjects = Integer.parseInt(tokens[0]);
			numTimeSteps = Integer.parseInt(tokens[1]);

			//Read next line and initialise row counter
			line = br.readLine();
			timeC = 0;
			objC = 0;

			while (line != null) {

				parseLine(line);
				line = br.readLine();
			}
			if (timeC != numTimeSteps) {
				System.out.print("Invalid environment file\n");
				System.exit(1);
			}
			problemLoaded = true;
			br.close();
		} catch (IOException e) {
			System.out.print(e);
		}

	}

	private void parseLine(String line){
		String delim = " ";
		String[] tokens = line.split(delim);
		int tokenCount = tokens.length;
		if (tokenCount == 2){
			timeC += 1;
		} else if (tokenCount == 7) {
			ObjConfig config = new ObjConfig();
			config.xpos = Double.parseDouble(tokens[3]);
			config.ypos = Double.parseDouble(tokens[4]);
			config.xvel = Double.parseDouble(tokens[5]);
			config.yvel = Double.parseDouble(tokens[6]);
			if (objC == 0){
				List<ObjConfig> configs = new ArrayList<ObjConfig>();
				configs.add(config);
				timesteps.put(timeC, configs);
			} else {
				List<ObjConfig> configs = timesteps.get(timeC);
				configs.add(config);
				timesteps.put(timeC, configs);
			}


			objC += 1;
			if (objC >= numObjects){
				objC = 0;
			}
		} else {
			System.out.print("Invalid environment file\n");
			System.exit(1);
		}
	}

	/**
	 * Returns whether a problem is currently loaded.
	 *
	 * @return whether a problem is currently loaded.
	 */
	public boolean problemLoaded() {
		return problemLoaded;
	}
}
