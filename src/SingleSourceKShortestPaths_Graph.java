import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Vector;

/*
 * this file implements the network (graph) and the algorithm.
 */

public class SingleSourceKShortestPaths_Graph {
	
	static final double Epsilon = 0.15f;  //the minimum weight of edges
	static final int Max_Nodes_In_A_Path = 23;  //equivalent to "h" in our paper  
		
	static double MAX_DISTANCE =  Double.MAX_VALUE ;  //default distance	
	
	int numGenes;//number of nodes in this graph
	int numEdges;//number of edges in this graph
	int tgnum;  //the target gene's number
	int cgnum;  //the causal gene's number
	Set<Integer> candidateCausalGenes = new HashSet<Integer>(); //the candidate genes' numbers
	
	Vector<TreeNode>[] directEdges;
	Importance[] importances;
	
	protected final static int numReadings = 5;
	
	SingleSourceKShortestPaths_Graph(){
	}
	
	SingleSourceKShortestPaths_Graph(int numGenes, int tg, int cg ){
		if(SingleSourceKShortestPaths.InferringTarget){
			//the program always starts at the node with the number tgnum, 
			//so in unknownTarget, we swap the target gene and causal gene 
			this.tgnum = cg;
			this.cgnum = tg;
		}
		else{
			this.tgnum = tg;
			this.cgnum = cg;
		}
		this.numGenes = numGenes;
		directEdges = new Vector[numGenes];
		importances = new Importance[numGenes];
		for(int i=0;i<numGenes;i++){
			directEdges[i] = new Vector<TreeNode>();
			importances[i] = new Importance(i);
		}	
	}

	//simulate the tree node and the path connecting to it.
	class TreeNode{
		int graphNodeNum; //the corresponding graph node number.
		double distance;
		int count; //the number of time the path appears in the pseudo tree.  this is only for K diverse short path algorithm
		
		TreeNode(int d, double distance){
			this.graphNodeNum = d;
			this.distance = distance;
			count = 0;
		}
		public boolean equals(Object o){
			TreeNode di = (TreeNode)o;
			return this.graphNodeNum == di.graphNodeNum;
		}
		public int hashCode(){  //used for equals
			return this.graphNodeNum;  
		}

	}
	
	//used to store all results for each node, including the path, the path distance, and the importance value 
	class Importance implements Comparable<Importance>{  //the importance value.
		int node; //the node number
		
		double distances[] = new double[SingleSourceKShortestPaths.KShortestPaths];  //sorted array, used to store top K distance 
		Vector<TreeNode> paths[] = new Vector[SingleSourceKShortestPaths.KShortestPaths];  //the corresponding paths of of distances[]
		double importance = 0; //the importance value

		
		Importance(){
		}
		
		Importance(int node){
			this.node = node;
			for(int i=0; i<SingleSourceKShortestPaths.KShortestPaths; i++){
				distances[i] = MAX_DISTANCE;
				paths[i] = new Vector<TreeNode>();
			}
		}

		boolean addPath(DistanceWalk pf){  //store the distance and the path
			for(int i=0; i<SingleSourceKShortestPaths.KShortestPaths; i++)
				if(pf.currentDistance < distances[i]){
					for(int j=SingleSourceKShortestPaths.KShortestPaths-1; j>i ; j--){
						distances[j] = distances[j-1];
						paths[j] = paths[j-1];
					}
					distances[i] = pf.currentDistance;
					paths[i] = new Vector<TreeNode>(pf.visitedNodes);
					return true;
				}
			return false;	
		}
		
		//calculate the importance value according to the paths' length
		double getImportance(){
			this.importance = 0;
			int pathIdx=0;

			for(; pathIdx < SingleSourceKShortestPaths.KShortestPaths ; pathIdx++){
				this.importance += 1/distances[pathIdx];
			}
			return (double)this.importance;
		}

		public int compareTo(Importance im){ //larger is better 
			if(this.importance > im.importance)
				return -1;
			else if(this.importance < im.importance)
				return 1;
			else
				return 0;
		}
		
	}
	
	//used to simulate the tree path; the purpose of implementing comparable is for the priority queue 
	class DistanceWalk implements Comparable<DistanceWalk> {
		Vector<TreeNode> visitedNodes;
		double currentDistance;		
		int currentNode;
		
		DistanceWalk(int startNode){
			visitedNodes = new Vector<TreeNode>();
			visitedNodes.add(new TreeNode(startNode,1));
			currentDistance = 0;
			currentNode = startNode;
		}
		DistanceWalk(DistanceWalk pf) {  //duplicate
			visitedNodes = new Vector<TreeNode>(pf.visitedNodes);
			currentDistance = pf.currentDistance;
			currentNode = pf.currentNode;
		}

		boolean goNext(TreeNode di) {  //return whether this tree path reach candidate causal genes
			
			//if(tgnum == di.graphNodeNum)
			//	System.out.println("error");
			currentDistance = currentDistance + di.distance;

			currentNode = di.graphNodeNum;
			visitedNodes.add(di);
			if(candidateCausalGenes.contains(di.graphNodeNum)) {
				return true;
			}
			return false;
		}
		
		public int compareTo(DistanceWalk dw) { //smaller is picked first 
			if(this.currentDistance < dw.currentDistance)
				return -1;
			else if(this.currentDistance > dw.currentDistance)
				return 1;
			else
				return 0;
		}

	}
	
	Vector<TreeNode> unvisitedNodes(DistanceWalk pf){
		Vector<TreeNode> nodes = new Vector<TreeNode>(this.directEdges[pf.currentNode]);  //shallow copy

		nodes.removeAll(pf.visitedNodes);
		return nodes;
	}


	
	void addDirectEdge(int source, int graphNodeNum, double weight){
		//if(weight > 1)
			//System.out.println("error: weight > 1");
		double distance = (double) ( -1*Math.log(weight)/Math.log(10) + 1);   // //convert to distance  +1 because each edge should has a minimum distance

		TreeNode d = new TreeNode(graphNodeNum,distance);
		this.directEdges[source].add(d);
		numEdges++;
	}

	void deleteDirectEdge(int source, int graphNodeNum){
		directEdges[source].remove(graphNodeNum);
		numEdges--;
	}

	void setCandidateCausalGenes(Set<Integer> candidates){
		this.candidateCausalGenes = candidates;
	}


	/*
	 * The main process called by the main function of SingleSourceKShortestPaths.
	 * We do not explicitly construct the pseudo tree. Instead, we store a series of flows beginning from the given node to simulate the pseudo tree,
	 * and store the paths for each node.
	 */
	ExecuteInfo flowStart() {
		
		long maxHeap = ExecuteInfo.memoryUsed();
		
		int source = this.tgnum;   
		long timeStart = System.currentTimeMillis();
		long timeToReadMem = 0;
		DistanceWalk firstTreePath = new DistanceWalk(source);

		PriorityQueue<DistanceWalk> treePaths = new PriorityQueue<DistanceWalk>();  //used to pick the candidate path with minimal distance.
		treePaths.add(firstTreePath);
		for(int i=1; i<directEdges[source].size(); i++){
			DistanceWalk newTreePath = new DistanceWalk(firstTreePath);
			if(!newTreePath.goNext(directEdges[source].elementAt(i)) || !SingleSourceKShortestPaths.StopCandidate) {
				treePaths.add(newTreePath);
			}
			storePath(newTreePath);	
		}
		//System.out.println("treePaths.size(): "+treePaths.size()+" "+directEdges[source].size());
		if(directEdges[source].size()>0){
			if(firstTreePath.goNext(directEdges[source].elementAt(0)) && SingleSourceKShortestPaths.StopCandidate)
				treePaths.remove(firstTreePath);
			storePath(firstTreePath);
		}
		
		int iterations = 0;
		int modValue = Math.round( ((float) numGenes) / ((float) numReadings) );
		while(!treePaths.isEmpty()){
			DistanceWalk treePath = treePaths.peek();				
			if( treePath.visitedNodes.size() >= Max_Nodes_In_A_Path){  //flow should be terminated based on thresholds, stop spreading
				treePaths.remove(treePath);
				continue;
			}
			Vector<TreeNode> destinations = unvisitedNodes(treePath);  //destinations are the unvisited nodes
			treePaths.remove(treePath);
			
			for(int i=0; i< destinations.size(); i++){
				DistanceWalk newTreePath = new DistanceWalk(treePath);
				boolean reachCausal = newTreePath.goNext(destinations.elementAt(i));
				if(storePath(newTreePath) ){   //storePath(newTreePath) is for optimization!!
					if((!reachCausal || !SingleSourceKShortestPaths.StopCandidate))						
						treePaths.add(newTreePath);     
				}
			}
			
			iterations++;
			if(iterations % modValue == 0) {
				long x = System.currentTimeMillis();
				maxHeap = Math.max(maxHeap, ExecuteInfo.memoryUsed());
				long y = System.currentTimeMillis();
				timeToReadMem += y-x;
			}

			
			//System.out.println("pf: "+pf.currentNode+" "+pf.currentDistance+" "+pf.visitedNodes.size() + " "+treePaths.size() + ";");
		}
		for(int i=0; i < numGenes  ; i++){   
			this.importances[i].getImportance();
		}		
		//Verify that testdatagen worked
//		System.out.println("n="+numGenes);
//		System.out.println("m="+numEdges);
		
		ExecuteInfo result = new ExecuteInfo((System.currentTimeMillis()-timeStart-timeToReadMem)/1000F, maxHeap);
		//result.timeToReadMem = timeToReadMem;
		return result;
		
	}


	boolean storePath(DistanceWalk pf){
		return this.importances[pf.currentNode].addPath(pf);
	}

	//output true only if it's in candidateCausal and the predict causal gene is correct
	boolean outputResult(PrintStream output, PrintStream ResultRankOutputStream){
		//sort importance value
		int rank = numGenes;
		Arrays.sort(this.importances);
		int numTotalGene = SingleSourceKShortestPaths.NumTotalGene;
		
		int highestRankCandidate = 999999, highestRankAmongCandidate = 999999;
		for(int i=0; i< numTotalGene ; i++){
			if( i < highestRankAmongCandidate && candidateCausalGenes.contains(this.importances[i].node)){
				highestRankAmongCandidate = i; 
				highestRankCandidate = this.importances[i].node;
			}
			if(this.importances[i].node != cgnum){
				output.print(/*"rank:"+i+" "+*/this.importances[i].node + " "+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].node)+" "+importances[i].importance);
				output.println();
				if(SingleSourceKShortestPaths.RESULT_SHOW_PATHS)
					for(int pathth = 0; pathth < SingleSourceKShortestPaths.KShortestPaths; pathth++){
						output.print((pathth+1)+"th path(distance:" + this.importances[i].distances[pathth]+"):");
						for(int pathNode=0; pathNode < this.importances[i].paths[pathth].size()-1;pathNode++) {
							output.print(this.importances[i].paths[pathth].elementAt(pathNode).graphNodeNum+"("+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].paths[pathth].elementAt(pathNode).graphNodeNum)+")"+">");
						}
						output.print(this.importances[i].node+"("+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].node)+")	");
						output.println();
					}
			}
			else {
				output.print(/*"rank:"+i+" "+*/this.importances[i].node + " "+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].node)+" "+importances[i].importance+"	*");
				rank = i;
				output.println();
				if(SingleSourceKShortestPaths.RESULT_SHOW_PATHS)
					for(int pathth = 0; pathth < SingleSourceKShortestPaths.KShortestPaths; pathth++){
						output.print((pathth+1)+"th path(distance:"+this.importances[i].distances[pathth]+"):");
						for(int pathNode=0; pathNode<this.importances[i].paths[pathth].size()-1;pathNode++) {
							output.print(this.importances[i].paths[pathth].elementAt(pathNode).graphNodeNum+"("+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].paths[pathth].elementAt(pathNode).graphNodeNum)+")"+">");
						}
						output.print(this.importances[i].node+"("+SingleSourceKShortestPaths.GeneInNetwork.get(this.importances[i].node)+")	");
						output.println();
					}
				
				
			}
		}
		if(highestRankAmongCandidate == rank){  //correct in CandidateCausal 
			ResultRankOutputStream.println(SingleSourceKShortestPaths.GeneInNetwork.get(tgnum)+" "+SingleSourceKShortestPaths.GeneInNetwork.get(cgnum)+" "+(rank+1)+" *");
			return true;
		}
		else {
			ResultRankOutputStream.println(SingleSourceKShortestPaths.GeneInNetwork.get(tgnum)+" "+SingleSourceKShortestPaths.GeneInNetwork.get(cgnum)+" "+(rank+1)+" "+highestRankCandidate+" "+SingleSourceKShortestPaths.GeneInNetwork.get(highestRankCandidate));
		}
		
		return false;
	}
}
