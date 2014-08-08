import java.io.PrintStream;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Vector;


public class SingleSourceKDiverseShortPaths_Graph extends SingleSourceKShortestPaths_Graph{
	int numNodesFoundKDivPath;
	ImportanceMultipleRun[] importances;  //the importance of each node
	
	static double MAX_DISTANCE =  Double.MAX_VALUE ;  //default distance	

	
	SingleSourceKDiverseShortPaths_Graph(int numGenes, int tg, int cg){
		super(numGenes, tg, cg);
		importances = new ImportanceMultipleRun[numGenes];
		for(int i=0;i<numGenes;i++){
			importances[i] = new ImportanceMultipleRun(i);
		}	
		numNodesFoundKDivPath = 0;
	}
	


	ExecuteInfo multipleRuns(){  //the process which is called by the main function of SingleSourceKDiverseShortPath
		long timeStart = System.currentTimeMillis();
		int numItr = 0;
		long maxHeap = 0;
		long timeToReadMem = 0;
		
		Random aRandom = new Random();

		while(numNodesFoundKDivPath < numGenes){  
			int numDeletedEdges = 0;
			ExecuteInfo stats = this.flowStart();  //start K shortest paths algorithm
			maxHeap = Math.max(stats.maxHeap, maxHeap);
			timeToReadMem += stats.timeToReadMem;
			
			for(int i=0; i < numGenes ; i++){
				importances[i].initialize(); 
			}
			//delete edges according to the counts of edges
			for(int i=0; i<numGenes; i++)
				for(int j=0; j < directEdges[i].size(); j++){
					int count = directEdges[i].elementAt(j).count;
					if(count > SingleSourceKDiverseShortPaths.KShortestPaths / 2){  //only count / SingleSourceKDiverseShortPaths.KShortestPaths > 1/2 would be deleted
						float prob = (float)count / SingleSourceKDiverseShortPaths.KShortestPaths;
						if(aRandom.nextFloat() < prob){
							this.deleteDirectEdge(i,j);
							j--;
							numDeletedEdges++;
						}
					}
				}
			//System.out.println("number of delete edges:"+numDeletedEdges);
			//the count of each node is initialized to 0;
			for(int i=0; i<numGenes; i++)
				for(int j=0; j < directEdges[i].size(); j++)
					directEdges[i].elementAt(j).count = 0;
			
			numItr++;
			if( numDeletedEdges < (double)numEdges/numGenes)  //if there are no deleted edges in this iteration, convergence and so stop the iteration
				break;

		}
//		System.out.println("total iterations:"+numItr);
		ExecuteInfo result = new ExecuteInfo((System.currentTimeMillis()-timeStart-timeToReadMem)/1000F, numItr);
		result.maxHeap = maxHeap;
		return result;
	}


	ExecuteInfo flowStart() {
		int source = this.tgnum;   
		long timeStart = System.currentTimeMillis();
		ExecuteInfo result = new ExecuteInfo();
		result.timeToReadMem = 0;
		result.maxHeap = ExecuteInfo.memoryUsed();

		PriorityQueue<DistanceWalk> pfs = new PriorityQueue<DistanceWalk>();
		DistanceWalk firstpf = new DistanceWalk(source);
		pfs.add(firstpf);
		for(int i=1; i<directEdges[source].size(); i++){
			DistanceWalk newpf = new DistanceWalk(firstpf);
			if(!newpf.goNext(directEdges[source].elementAt(i)) || !SingleSourceKDiverseShortPaths.StopCandidate)
				pfs.add(newpf);
			storePath(newpf);	
		}
		if(directEdges[source].size() > 0){
			if(firstpf.goNext(directEdges[source].elementAt(0)) && SingleSourceKDiverseShortPaths.StopCandidate)
				pfs.remove(firstpf);
			storePath(firstpf);
		}
		
		int iterations = 0;
		int modValue = Math.round( ((float) numGenes) / ((float) numReadings) );
		
		while(!pfs.isEmpty()){
			DistanceWalk pf = pfs.peek();				
			if( pf.visitedNodes.size() >= Max_Nodes_In_A_Path){  //flow should be terminated based on thresholds, stop spreading
				pfs.remove(pf);
				continue;
			}
			Vector<TreeNode> destinations = unvisitedNodes(pf);
			pfs.remove(pf);

			for(int i=0; i< destinations.size(); i++){
				DistanceWalk newpf = new DistanceWalk(pf);
				boolean reachCausal = newpf.goNext(destinations.elementAt(i));
				if(storePath(newpf) ){   //storePath(newpf) is for optimization!!
					destinations.elementAt(i).count++; 
					if((!reachCausal || !SingleSourceKDiverseShortPaths.StopCandidate))						
						pfs.add(newpf);     
				}
			}
			
			iterations++;
			if(iterations % modValue == 0) {
				long x = System.currentTimeMillis();
				result.maxHeap = Math.max(result.maxHeap, ExecuteInfo.memoryUsed());
				long y = System.currentTimeMillis();
				result.timeToReadMem += y-x;
			}
		}
		//System.out.println("total exe. time (secs):" + (System.currentTimeMillis()-timeStart)/1000F);	
		for(int i=0; i < numGenes  ; i++){   
			this.importances[i].getImportance();
		}
		return result;
		
	}

	boolean storePath(DistanceWalk pf){
		return this.importances[pf.currentNode].addPath(pf);
	}


	//store all results for each node, including the path, the path distance, and the importance value 
	class ImportanceMultipleRun extends Importance{

		int numDiversePaths;  //this is a object variable instead of a local variable in class Graph 
		double totalDistance; //total distance of K diverse Paths
		double diversity;

		Vector<TreeNode> diversePaths[] = new Vector[SingleSourceKDiverseShortPaths.KShortestPaths];  //the diverse paths
		double[] diverseDistance = new double[SingleSourceKDiverseShortPaths.KShortestPaths];

		ImportanceMultipleRun(int node){
			super(node);
			numDiversePaths = 0;
			totalDistance = 0;
			this.importance = 0;
			for(int i=0; i<SingleSourceKDiverseShortPaths.KShortestPaths; i++){
				diversePaths[i] = new Vector<TreeNode>();
				diverseDistance[i] = MAX_DISTANCE;
			}
		}

		void initialize(){// call it each run   delete all distance and paths information from last runs
			distances = new double[SingleSourceKDiverseShortPaths.KShortestPaths];  
			paths = new Vector[SingleSourceKDiverseShortPaths.KShortestPaths];
			for(int i=0; i<SingleSourceKDiverseShortPaths.KShortestPaths; i++){
				distances[i] = MAX_DISTANCE;
				paths[i] = new Vector<TreeNode>();
			}
		}

		float getDiversity(){	//output the diversity of K diverse paths after finding them.
			return 0; //total distinct edges in K diverse paths / sum(# edges) of all paths
		}

		//calculate the importance value.
		double getImportance(){
			if(numDiversePaths == SingleSourceKDiverseShortPaths.KShortestPaths)  //enough diverse paths are found, do not need to add more path and update improtance value
				return (double)this.importance;

			double diverse_paths = 0;
			int pathIdx=0;
			for(; pathIdx < SingleSourceKDiverseShortPaths.KShortestPaths ; pathIdx++){
				
				//calculate the diversity = # of edges not appear in previous paths / total edgesin this path
				if(paths[pathIdx].size() == 0){  //only happen for the start node or there are not enough K paths
					for(int i=numDiversePaths; i < SingleSourceKDiverseShortPaths.KShortestPaths; i++)
						diversePaths[i] = paths[pathIdx];
					totalDistance += MAX_DISTANCE * (SingleSourceKDiverseShortPaths.KShortestPaths-numDiversePaths);
					numDiversePaths = SingleSourceKDiverseShortPaths.KShortestPaths;  //won't calculate importance in the future
					numNodesFoundKDivPath++; //System.out.println("paths[pathIdx].size() == 0");
					break;
				}
				
				boolean[] distinctEdges = new boolean[paths[pathIdx].size()-1];
				for(int i=0; i<distinctEdges.length ; i++)
					distinctEdges[i] = true;
				for(int edgeIdx = 0; edgeIdx < paths[pathIdx].size()-1; edgeIdx ++){
					for(int divPathIdx = 0; divPathIdx < numDiversePaths; divPathIdx++){
						for(int edgeIdx2 = 0; edgeIdx2 < diversePaths[divPathIdx].size()-1; edgeIdx2 ++){
							if(paths[pathIdx].elementAt(edgeIdx).graphNodeNum == diversePaths[divPathIdx].elementAt(edgeIdx2).graphNodeNum && paths[pathIdx].elementAt(edgeIdx+1).graphNodeNum == diversePaths[divPathIdx].elementAt(edgeIdx2+1).graphNodeNum){  //the same edge is found in 1~(pathIdx-1) shortest diverse paths
								distinctEdges[edgeIdx] = false;	
								break;
							}
						}
						if(!distinctEdges[edgeIdx])
							break;	
					}
				}
				int tempCount = 0;
				for(int i=0; i<distinctEdges.length ; i++)
					if(distinctEdges[i])
						tempCount++;
				double diversity = (double)tempCount / distinctEdges.length;  //this is the diversity of a specific path
				//System.out.println("Div:"+diversity+" ");
				
				//if the path's diversity exceeds the threshold, store the path and update the importance value 
				if(diversity >= SingleSourceKDiverseShortPaths.DiversityThreshold){
					diversePaths[numDiversePaths] = paths[pathIdx];
					diverseDistance[numDiversePaths] = distances[pathIdx];
					numDiversePaths++;
					this.importance += 1/distances[pathIdx];  //only count the paths are considered as diverse paths
					this.totalDistance += distances[pathIdx];				

					if(numDiversePaths == SingleSourceKDiverseShortPaths.KShortestPaths){  // do not need to calculate following paths, since the number
						numNodesFoundKDivPath++;
						break;
					}
					
				}
				
				
			}
			//System.out.print(pathIdx+" "+diverse_paths+";");
			return (double)this.importance;
		}

	}

	boolean outputResult(PrintStream output, PrintStream ResultRankOutputStream){

		//output average total distance versus average diversity of k paths
		double totalDiversity = 0;
		double totalDistance = 0;
		int connectedGenes = 0;
		double totalImps = 0;
		for(int nodeIdx=0; nodeIdx < numGenes ; nodeIdx++){
			if(nodeIdx != this.tgnum){
				
				//diversity of K paths from the source node to i node = total distinct edges / sum(total edges) of all paths
				int numDistinctEdges = Math.max(0, importances[nodeIdx].diversePaths[0].size()-1); //.size() is the number of nodes, so .size()-1 is the number of edges
				int totalEdges = Math.max(0, importances[nodeIdx].diversePaths[0].size()-1);
				for(int pathIdx = 1; pathIdx < SingleSourceKDiverseShortPaths.KShortestPaths; pathIdx++ ){
					totalEdges += Math.max(0, importances[nodeIdx].diversePaths[pathIdx].size()-1);
					for(int edgeIdx = 0; edgeIdx < importances[nodeIdx].diversePaths[pathIdx].size()-1; edgeIdx ++){
						boolean distinct = true;
						for(int pathIdx2 = 0; pathIdx2 < pathIdx; pathIdx2++){
							for(int edgeIdx2 = 0; edgeIdx2 < importances[nodeIdx].diversePaths[pathIdx2].size()-1; edgeIdx2 ++){
								if(importances[nodeIdx].diversePaths[pathIdx].elementAt(edgeIdx).graphNodeNum == importances[nodeIdx].diversePaths[pathIdx2].elementAt(edgeIdx2).graphNodeNum && importances[nodeIdx].diversePaths[pathIdx].elementAt(edgeIdx+1).graphNodeNum == importances[nodeIdx].diversePaths[pathIdx2].elementAt(edgeIdx2+1).graphNodeNum) //the same edge
									distinct = false;
							}
						}
						if(distinct)
							numDistinctEdges++;
					}
				}
				if(totalEdges > 0){ //eliminate disconnected nodes	
					importances[nodeIdx].diversity = (double)numDistinctEdges/totalEdges;
					totalDiversity += importances[nodeIdx].diversity ;
					totalDistance += importances[nodeIdx].totalDistance ;
					totalImps += importances[nodeIdx].importance;
					
					connectedGenes++;
				}
			}
		}	
//		System.out.println("average importance value:"+ totalImps/connectedGenes+"  ;"+connectedGenes);
//		System.out.println("average diversity:"+ totalDiversity/connectedGenes );
		//System.out.println("average total distance:"+ totalDistance/connectedGenes );
		

		//sort importance value
		int rank = numGenes;
		Arrays.sort(this.importances);

		int highestRankCandidate = 999999, highestRankAmongCandidate = 999999;
		double maxImp = 0;
		for(int i=0; i< SingleSourceKDiverseShortPaths.NumTotalGene ; i++){
			if(i == 0)
				maxImp = importances[i].importance;
			if( i < highestRankAmongCandidate && candidateCausalGenes.contains(this.importances[i].node)){
				highestRankAmongCandidate = i; 
				highestRankCandidate = this.importances[i].node;
			}
			if(this.importances[i].node != cgnum){
				double impRatio = importances[i].importance/maxImp;
				output.print(/*"rank:"+i+" "+*/this.importances[i].node + " "+SingleSourceKDiverseShortPaths.GeneInNetwork .get(this.importances[i].node)+" "+importances[i].importance+" div:"+importances[i].diversity	  );
				output.println();
				if(SingleSourceKDiverseShortPaths.RESULT_SHOW_PATHS)
					for(int pathth = 0; pathth < SingleSourceKDiverseShortPaths.KShortestPaths; pathth++){
						output.print((pathth+1)+"th path(distance:" + this.importances[i].diverseDistance[pathth]+"):");
						for(int pathNode=0; pathNode < this.importances[i].diversePaths[pathth].size()-1;pathNode++)
							output.print(this.importances[i].diversePaths[pathth].elementAt(pathNode).graphNodeNum+"("+SingleSourceKDiverseShortPaths.GeneInNetwork.get(this.importances[i].diversePaths[pathth].elementAt(pathNode).graphNodeNum)+")"+">");
						output.print(this.importances[i].node+"("+SingleSourceKDiverseShortPaths.GeneInNetwork.get(this.importances[i].node)+")	");
						output.println();
					}
			}
			else{
				output.print(/*"rank:"+i+" "+*/this.importances[i].node + " "+SingleSourceKDiverseShortPaths.GeneInNetwork.get(this.importances[i].node)+" "+importances[i].importance+"	*");
				rank = i;
				output.println();
				if(SingleSourceKDiverseShortPaths.RESULT_SHOW_PATHS)
					for(int pathth = 0; pathth < SingleSourceKDiverseShortPaths.KShortestPaths; pathth++){
						output.print((pathth+1)+"th path(distance:"+this.importances[i].diverseDistance[pathth]+"):");
						for(int pathNode=0; pathNode<this.importances[i].diversePaths[pathth].size()-1;pathNode++)
							output.print(this.importances[i].diversePaths[pathth].elementAt(pathNode).graphNodeNum+"("+SingleSourceKDiverseShortPaths.GeneInNetwork.get(this.importances[i].diversePaths[pathth].elementAt(pathNode).graphNodeNum)+")"+">");
						output.print(this.importances[i].node+"("+SingleSourceKDiverseShortPaths.GeneInNetwork.get(this.importances[i].node)+")	");
						output.println();
					}
				
				
			}
		}
		if(highestRankAmongCandidate == rank){  //correct
			ResultRankOutputStream.println(SingleSourceKDiverseShortPaths.GeneInNetwork.get(tgnum)+" "+SingleSourceKDiverseShortPaths.GeneInNetwork.get(cgnum)+" "+(rank+1)+" *");
			return true;
		}
		else
			ResultRankOutputStream.println(SingleSourceKDiverseShortPaths.GeneInNetwork.get(tgnum)+" "+SingleSourceKDiverseShortPaths.GeneInNetwork.get(cgnum)+" "+(rank+1)+" "+highestRankCandidate+" "+SingleSourceKDiverseShortPaths.GeneInNetwork.get(highestRankCandidate));

		return false;
	}
}
