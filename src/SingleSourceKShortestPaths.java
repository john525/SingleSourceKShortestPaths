import java.awt.Graphics;
import java.awt.Image;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;

/*
 * This program implements the single-source K shortest paths algorithm. More detail can be found in: Yu-Keng Shih and Srinivasan Parthasarathy, A single-source k shortest paths algorithm to infer regulatory pathways in a gene network, to appear in ISMB 2012
 * args[0] is the file indicating all gold standards. The format is described in readme.
 * args[1] is the gene number and gene name mapping. The format is described in readme.
 * args[2] is K, the number of shortest paths
 * args[3] is the choice of mode (1 for unknownCausal, 2 for unknownTarget, 3 for candidateCausal)
 * output is a set of weighted graph for each target gene.
 * args[4] is the GraphDir = InvGraphDir.

 * The program is made by Yu-Keng Shih from the Ohio State University
 * Last update: Aug 5, 2012

 * If you have any question or find any bug, please contact Yu-Keng Shih: shihy@cse.ohio-state.edu

 * 
 */

public class SingleSourceKShortestPaths {
	
	static String l = System.getProperty("file.separator");

	//modify the directory if necessary
	String GraphDir = ".."+l+"sample_input_graphs"+l;  //the directory containing graphs for unknownCausal and candidateCausal
	String InvGraphDir = ".."+l+"sample_input_graphs"+l;  //the directory containing graphs for unknownTarget
	
	String ResultDirectory;
	static final String RESULT_DIR_UNKNOWNCAUSAL = "result_causal"+l;
	static final String RESULT_DIR_UNKNOWNTARGET = "result_target"+l;
	static final String RESULT_DIR_CANDIDATECAUSAL = "result_candidate"+l;
	
	
	static int KShortestPaths;    //equivalent to "k": the number of paths
	static final boolean RESULT_SHOW_PATHS = false;    //the option for output. If true, it will output the path detail; otherwise, it only output the importance value.
	
	static int NumTotalGene = 0; //the number of total genes     
	
	static int NumCandidateCausalGene = 10;
	
	static boolean StopCandidate;   //the tree path will stop expand when the tree path arrive a candidate causal gene. It's true only in CandidateCausal
	static boolean InferringTarget;	//false for inferring causal gene
	
	static Hashtable<Integer,String> nodeNameNetwork = new Hashtable<Integer,String>();	//map the number for gene to the gene name
	
	static long heapStart;
	
	public static void main(String[] args) {
		doTesting(new SingleSourceKShortestPaths());
	}
	
	public static void doTesting(SingleSourceKShortestPaths algo) {
		JFrame frame = new JFrame("K Shortest Paths Algorithm");
		frame.setSize(669,470);
		frame.setLocationRelativeTo(null);
		
		JTextArea text = new JTextArea();
		text.setText("Loading...");
		frame.add(text);
		frame.setVisible(true);
		
		final int numTrials = 3;
		double[][][] memData = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length][numTrials];
		double[][][] timeData = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length][numTrials];

		TestDataGenerator gen = new TestDataGenerator();
		
		for(int i=0; i<numTrials; i++) {
			gen.go();
			
			for(int size=0; size<TestDataGenerator.sizes.length; size++) {
				for(int deg=0; deg<TestDataGenerator.avgDegrees.length; deg++) {
					String loc = "test_data";
					String subFolder = loc + l + "n=" + TestDataGenerator.sizes[size] + ", deg=" + TestDataGenerator.avgDegrees[deg];
					String[] x = {subFolder+l+"allgraphs_all", subFolder+l+"allProteins", "5", "2", subFolder+l, "0.5"};
					System.out.println("Attempting n=" + TestDataGenerator.sizes[size] + " and deg=" + TestDataGenerator.avgDegrees[deg]);
					
					
					try {
						ExecuteInfo res = algo.runAlgorithm( x );
						memData[size][deg][i] = res.maxHeap;
						timeData[size][deg][i] = res.time;
						System.out.println(res);
						System.out.println();
					}
					catch(Exception e) {//IO or maybe NullPointer
						e.printStackTrace();
					}
				}
			}
		}
		
		System.out.println("*****************************************");
		System.out.println("**********RESULTS START HERE*************");
		System.out.println("*****************************************");
		double[][] memAverages = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length];
		double[][] memSD = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length];
		double[][] timeAverages = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length];
		double[][] timeSD = new double[TestDataGenerator.sizes.length][TestDataGenerator.avgDegrees.length];

		
		StringBuilder memBuilder = new StringBuilder(), timeBuilder = new StringBuilder();
		for(int size=0; size<TestDataGenerator.sizes.length; size++) {
			for(int deg=0; deg<TestDataGenerator.avgDegrees.length; deg++) {
				for(int i=0; i<numTrials; i++) {
					memAverages[size][deg] += memData[size][deg][i];
					timeAverages[size][deg] += timeData[size][deg][i];
				}
				memAverages[size][deg] /= (double)numTrials;
				timeAverages[size][deg] /= (double)numTrials;

				for(int i=0; i<numTrials; i++) {
					memSD[size][deg] += Math.pow(memAverages[size][deg] - memData[size][deg][i], 2);
					timeSD[size][deg] += Math.pow(timeAverages[size][deg] - timeData[size][deg][i], 2);
				}
				memSD[size][deg] /= (double)numTrials;
				memSD[size][deg] = Math.sqrt(memSD[size][deg]);
				
				timeSD[size][deg] /= (double)numTrials;
				timeSD[size][deg] = Math.sqrt(timeSD[size][deg]);

				memBuilder.append("n="+TestDataGenerator.sizes[size]+", deg="+TestDataGenerator.avgDegrees[deg]+": "+memAverages[size][deg]+"�"+memSD[size][deg]+"\n");
				timeBuilder.append("n="+TestDataGenerator.sizes[size]+", deg="+TestDataGenerator.avgDegrees[deg]+": "+timeAverages[size][deg]+"�"+timeSD[size][deg]+"\n");
			}
		}

		System.out.println(memBuilder);
		System.out.println();
		System.out.println(timeBuilder);
		try {
			long timeStamp = System.currentTimeMillis();
			File memFile = new File("k_memory_"+timeStamp);
			if(!memFile.exists()) memFile.createNewFile();
			File timeFile = new File("k_time_"+timeStamp);
			if(!timeFile.exists()) timeFile.createNewFile();
			PrintStream stream = new PrintStream(memFile);
			stream.println(memBuilder.toString());
			stream.close();
			PrintStream stream2 = new PrintStream(timeFile);
			stream2.println(timeBuilder.toString());
			stream2.close();
			text.setText("TIME:\n" + timeBuilder.toString() + "\nMEMORY:\n" + memBuilder.toString());
			frame.repaint();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public ExecuteInfo runAlgorithm(String args[]) throws IOException {		
		heapStart = ExecuteInfo.memoryUsed();
		
		NumTotalGene = 0;
		nodeNameNetwork = new Hashtable<Integer,String>();

		KShortestPaths = Integer.valueOf(args[2]);
		GraphDir = InvGraphDir = args[4];
		int mode = Integer.valueOf(args[3]);
		
		if(mode == 1){  //unknownCausal
			StopCandidate = false;
			InferringTarget = false;
			ResultDirectory = RESULT_DIR_UNKNOWNCAUSAL;
			NumCandidateCausalGene = 1;
		}
		else if(mode == 2){   //unknownTarget
			StopCandidate = false;
			InferringTarget = true;
			ResultDirectory = RESULT_DIR_UNKNOWNTARGET;
			NumCandidateCausalGene = 1;
		}
		else{  //CandidateCausal
			StopCandidate = true;
			InferringTarget = false;
			ResultDirectory = RESULT_DIR_CANDIDATECAUSAL;
		}
		
		GraphDir = InvGraphDir = args[4];

		File genelistFile = new File(args[1]);
		BufferedReader genelistReader = new BufferedReader(new FileReader(genelistFile));  
		
		/*Read the gene name and gene number mapping*/
		while(genelistReader.ready()){
			String geneName = genelistReader.readLine();
			nodeNameNetwork.put( NumTotalGene , geneName);
			NumTotalGene ++;
		}
		genelistReader.close();


		File allgoldFile = new File(args[0]);
		BufferedReader allgoldReader = new BufferedReader(new FileReader(allgoldFile));
		allgoldReader.readLine(); //Target_gene Target_gene_num Causal_gene Causal_gene_num
		
		File resultRank = new File(ResultDirectory+"rank_K=" + KShortestPaths+"_shortestPaths");
		if(!resultRank.exists()) resultRank.createNewFile();
		PrintStream ResultRankOutputStream = new PrintStream(resultRank);

		int correctInferCausal = 0;	//only work in CandidateCausal
		int totalNumCausal = 0;
		double totaltime = 0;
		
		ExecuteInfo res = null;
		
		while(allgoldReader.ready()){  //read each gold standard and process it
			String aLine = allgoldReader.readLine();
			String[] infos = aLine.split(" ");
			String tg = infos[0];  //target gene
			String cg = infos[2];  //causal gene
			int tgNum = Integer.valueOf(infos[1]);
			int cgNum = Integer.valueOf(infos[3]);
						
			File aGraphFile = new File(GraphDir+tg+"_"+cg);
			if(InferringTarget){
				GraphDir = InvGraphDir;
				aGraphFile = new File(GraphDir+tg+"_"+cg+"_inv");
			}
			BufferedReader aGraphReader = new BufferedReader(new FileReader(aGraphFile));
//			System.out.println("start processing "+totalNumCausal+"th file:" +aGraphFile.getName() +" K="+ KShortestPaths+" mode:"+mode);
			
			//candidate causal genes
			Set<Integer> candidateCausalGenes = new HashSet<Integer>();
			String[] infos2 = aGraphReader.readLine().split(" ");
			for(int i=0; i<NumCandidateCausalGene && i<infos2.length ; i++) {
				if(!infos2[i+2].equals("null")) {
					candidateCausalGenes.add(Integer.valueOf(infos2[i+2]));
				}
			}
			
			SingleSourceKShortestPaths_Graph aGraph = new SingleSourceKShortestPaths_Graph(NumTotalGene,tgNum,cgNum);
			aGraph.setCandidateCausalGenes(candidateCausalGenes);
			while(aGraphReader.ready()){
				String aLine2 = aGraphReader.readLine();
				String aEdgeInfo[] = aLine2.split(" ");
				aGraph.addDirectEdge(Integer.valueOf(aEdgeInfo[0]), Integer.valueOf(aEdgeInfo[1]), Float.valueOf(aEdgeInfo[2]));	
			}
			
			//FIXME only works if there's just one gold standard
			res = aGraph.flowStart();  //main process is here!
			double thistime = res.time;
			res.maxHeap = res.maxHeap - heapStart;
			
//			System.out.println("exe. time:"+thistime+" sec");
			totaltime += thistime;
			
			//Testing output code starts
			PrintStream statsOutput = new PrintStream(GraphDir + "algo stats(" + aGraphFile.getName() + ", Mode " + mode + ")");
			statsOutput.println("Time: " + thistime + " secs");
			statsOutput.println("Max Memory: " + res.maxHeap/1000000L + "MB");
			//Testing output code ends
			
//			System.out.println("Time: " + thistime + " secs");
//			System.out.println("Max Memory: " + res.maxHeap/1000000L + "MB");
//			System.out.println();
			
			File f = new File(ResultDirectory, "K="+KShortestPaths+"_shortestPaths"+"_tg="+tg+"_cg="+cg);
			if(!f.exists()) f.createNewFile();
			PrintStream resultOutputStream = new PrintStream(f);
			if(aGraph.outputResult(resultOutputStream,ResultRankOutputStream))
				correctInferCausal++;
			totalNumCausal++;
			resultOutputStream.close();
			aGraphReader.close();
		}
//		System.out.println("accuracy:"+correctInferCausal+"/"+totalNumCausal+" = "+((double)correctInferCausal/(double)totalNumCausal));
//		System.out.println("avg. execution time (secs):" + totaltime/(double)totalNumCausal );
		ResultRankOutputStream.close();
		allgoldReader.close();
		
		return res;
	}
	


}
