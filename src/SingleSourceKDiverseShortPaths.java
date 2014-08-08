import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;


/*
 * This program implements the single-source K shortest paths algorithm. More detail can be found in: Yu-Keng Shih and Srinivasan Parthasarathy,
 * A single-source k shortest paths algorithm to infer regulatory pathways in a gene network, to appear in ISMB 2012
 * args[0] is the file indicating all gold standards. The format is described in readme.
 * args[1] is the gene number and gene name mapping. The format is described in readme.
 * args[2] is K, the number of shortest paths
 * args[3] is the choice of mode (1 for unknownCausal, 2 for unknownTarget, 3 for candidateCausal)
 * output is a set of weighted graph for each target gene.  
 * args[4] is GraphDir = InvGraphDir
 * args[5] is lambda, the diversity threshold to determine whether a path is diverse or not.
 * the code is similar to SingleSourceKShortestPaths. So please refer to the comment in SingleSourceKShortestPaths.java

 * The program is made by Yu-Keng Shih from the Ohio State University
 * Last update: Aug 5, 2012

 * If you have any question or find any bug, please contact Yu-Keng Shih: shihy@cse.ohio-state.edu

 * 
 */

public class SingleSourceKDiverseShortPaths extends SingleSourceKShortestPaths {

	static double DiversityThreshold; 	
	
	public static void main(String[] args) {
		doTesting(new SingleSourceKDiverseShortPaths());
	}
	
	@Override
	public ExecuteInfo runAlgorithm(String args[]) throws IOException {
		long initialHeap = ExecuteInfo.memoryUsed(), maxHeap = initialHeap;
		
		KShortestPaths = Integer.valueOf(args[2]);
		DiversityThreshold = Double.valueOf(args[5]);
		GraphDir = InvGraphDir = args[4];
		int mode = Integer.valueOf(args[3]);
		if(mode == 1){
			StopCandidate = false;
			InferringTarget = false;
			ResultDirectory = RESULT_DIR_UNKNOWNCAUSAL;
			NumCandidateCausalGene = 1;
		}
		else if(mode == 2){
			StopCandidate = false;
			InferringTarget = true;
			ResultDirectory = RESULT_DIR_UNKNOWNTARGET ;
			NumCandidateCausalGene = 1;
		}
		else{
			StopCandidate = true;
			InferringTarget = false;
			ResultDirectory = RESULT_DIR_CANDIDATECAUSAL;
		}


		File genelistFile = new File(args[1]);
		BufferedReader genelistReader = new BufferedReader(new FileReader(genelistFile));  
		
		/*Read all genes*/
		while(genelistReader.ready()){
			String geneName = genelistReader.readLine();
			GeneInNetwork.put( NumTotalGene , geneName);
			NumTotalGene ++;
		}
		genelistReader.close();

		File allgoldFile = new File(args[0]);
		BufferedReader allgoldReader = new BufferedReader(new FileReader(allgoldFile));
		allgoldReader.readLine(); //Target_gene Target_gene_num Causal_gene Causal_gene_num
		PrintStream ResultRankOutputStream = new PrintStream(ResultDirectory+"rank_divShortPaths"+"_K="+KShortestPaths+"_lambda="+DiversityThreshold);

		int correctInferCausal = 0;
		int totalNumCausal = 0;
		double totaltime = 0, timeToReadMem = 0;
		int totalItr = 0;

		while(allgoldReader.ready()){
			String aLine = allgoldReader.readLine();
			String[] infos = aLine.split(" ");
			String tg = infos[0];
			String cg = infos[2];
			int tgNum = Integer.valueOf(infos[1]);
			int cgNum = Integer.valueOf(infos[3]);				

			
			File aGraphFile = new File(GraphDir+tg+"_"+cg);
			if(InferringTarget){
				GraphDir = InvGraphDir;
				aGraphFile = new File(GraphDir+tg+"_"+cg+"_inv");
			}
			BufferedReader aGraphReader = new BufferedReader(new FileReader(aGraphFile));
			System.out.println("start processing "+totalNumCausal+"th file:" +aGraphFile.getName() +
					" K="+ KShortestPaths+" DiversityThreshold="+DiversityThreshold+" mode:"+mode+" div2Run threshold"+DiversityThreshold);

			Set<Integer> candidateCausalGenes = new HashSet<Integer>();
			String infos2[] = aGraphReader.readLine().split(" ");
			for(int i=0; i<NumCandidateCausalGene && i<infos2.length ; i++) {
				if(!infos2[i+2].equals("null")) {
					candidateCausalGenes.add(Integer.valueOf(infos2[i+2]));
				}
			}
			
			SingleSourceKDiverseShortPaths_Graph aGraph = new SingleSourceKDiverseShortPaths_Graph(NumTotalGene, tgNum, cgNum);
			aGraph.setCandidateCausalGenes(candidateCausalGenes);
			while(aGraphReader.ready()){
				String aLine2 = aGraphReader.readLine();
				String aEdgeInfo[] = aLine2.split(" ");
				aGraph.addDirectEdge(Integer.valueOf(aEdgeInfo[0]), Integer.valueOf(aEdgeInfo[1]), Float.valueOf(aEdgeInfo[2]));	
			}

			ExecuteInfo executeInfo = aGraph.multipleRuns();  //main process!
			System.out.println("exe. time:" + executeInfo.time);
			totaltime += executeInfo.time;
			totalItr += executeInfo.numItr;
			maxHeap = Math.max(executeInfo.maxHeap, maxHeap);
			//timeToReadMem += executeInfo.timeToReadMem;

			PrintStream resultOutputStream = new PrintStream(ResultDirectory+"K="+KShortestPaths+"_divShortPaths_lambda="+DiversityThreshold+"_tg:"+tg+"_cg:"+cg);
			if(aGraph.outputResult(resultOutputStream, ResultRankOutputStream)) {
				correctInferCausal++;
			}
			totalNumCausal++;
			resultOutputStream.close();
			
		}  
//		System.out.println("accuracy:"+correctInferCausal+"/"+totalNumCausal+" = "+((double)correctInferCausal/(double)totalNumCausal));
//		System.out.println("avg. execution time (secs):" + totaltime);
//		System.out.println("avg. number of iterations:" + totalItr/(double)totalNumCausal );
		ResultRankOutputStream.close();
		
		ExecuteInfo result = new ExecuteInfo(totaltime, totalItr);
		result.maxHeap = maxHeap - initialHeap;
		return result;
	}
}
