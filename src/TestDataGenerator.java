import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;


public class TestDataGenerator {
	
	public static final int[] sizes = {10000, 50000, 250000, 1000000};
	//public static final float[] densities = {0.001f, 0.002f, 0.005f};
	public static final int[] avgDegrees = {5, 10, 20};
	public static final Random rand = new Random();
		
	void randomize(File f, int size, int deg) throws IOException {
		//true = autoflush
		PrintStream output = new PrintStream(new FileOutputStream(f), true);
		
		String header = size + " 0 1";
		output.println(header);
//		long v = size;
//		System.out.println("v=" + v);
//		int m = (int)Math.round(density * (v*(v-1)) );
//		//System.out.println("m=" + m);
//		BloomFilter edgeSet = new BloomFilter(m);//Can make a hash table for increased accuracy
//		for(int i=0; i<m; i++) {
//			System.out.println("i="+i);
//			Edge edge;
//			do {
//				edge = new Edge(rand.nextInt(size), rand.nextInt(size));
//			} while(edgeSet.contains(edge) || edge.head == edge.tail);
//			edgeSet.add(edge);
//			String text = edge.head + " " + edge.tail + " " + (1.0F-rand.nextFloat());
//			output.println(text);
//			if( i % 1000 == 0) {
//				System.out.println(size+": "+i);
//			}
//		}
		
		//FIXME Average degree isn't set right.
		//Maybe have an array with length "size," keep making random edges until you have deg*size edges, or something like that
		//Or keep this code, and modify it by adding an array to keep track of edge sizes, so you don't add more after it already has enough
		//(Right now, average outdegree is set to deg, but you are forgetting that indegree is added on too)
		
//		byte[] nodeDegrees = new byte[size];
//		for(int i=0; i<nodeDegrees.length; i++) {
//			nodeDegrees[i] = 0;
//		}
//		
//		//FIXME boundary conditions at end
//		for(int node = 0; node < size-1; node++) {
//			while(nodeDegrees[node] < deg) {
//				
//				if(Math.random() < Math.exp(-Math.abs(deg-nodeDegrees[node]))) break;
//				
//				int otherNode;
//				do {
//					otherNode = node + 1 + rand.nextInt(size - (node+1));
//				} while (nodeDegrees[otherNode] >= deg);
//				
//				nodeDegrees[node]++;
//				nodeDegrees[otherNode]++;
//				output.println(node + " " + otherNode + " " + (1.0-Math.random()));
//			}
//		}
		
		ArrayList<Integer>[] edges = new ArrayList[size];
		for(int i=0; i<edges.length; i++) {
			edges[i] = new ArrayList<Integer>();
		}
		long totalNumEdges = 0;
		String x;
		while (((float)2*totalNumEdges)/((float)size) < deg) {
			int a = rand.nextInt(size);
			int b = rand.nextInt(size-1);
			if(b >= a) b++;
			if(edges[a].contains(b)) {
				continue;
			}
			else {
				edges[a].add(b);
				x = a + " " + b + " " + (1.0-Math.random());
				output.println(x);
				totalNumEdges++;
			}
		}
		
		output.close();
	}
	
	public static void main(String[] args) {
		new TestDataGenerator().go();
	}
	
	void deleteFolder(File f) {
		if(f.isDirectory()) {
			for(File subFile : f.listFiles()) {
				deleteFolder(subFile);
			}
			f.delete();
		}
		else {
			f.delete();
		}
	}
	
	public void go() {
		
		//heapA = ExecuteInfo.memoryUsed();
		
		File testFolder = new File("test_data");
		if(testFolder.exists()) {
			deleteFolder(testFolder);
		}
		testFolder.mkdir();
				
		int count = 0;
		
		for(int size: sizes) {
			for(int deg : avgDegrees) {
				File graphFolder = new File(testFolder, "n="+size+", deg="+deg);
				if(!graphFolder.exists()) {
					graphFolder.mkdir();
				}
				File geneList = new File(graphFolder, "allProteins");
				File graphList = new File(graphFolder, "allgraphs_all");
				File unknownCausal = new File(graphFolder, "G0_G1");
				File unknownTarget = new File(graphFolder, "G0_G1_inv");
				BufferedWriter writer;
				
				try {
					writer = new BufferedWriter(new FileWriter(geneList));
					for(int geneNum=0; geneNum<size; geneNum++) {
						writer.write("G"+geneNum+"\n");
					}
					writer.flush();
					writer.close();
					
					writer = new BufferedWriter(new FileWriter(graphList));
					writer.write("Target_gene Target_gene_num Causal_gene Causal_gene_num\n");
					writer.write("G0 0 G1 1");
					writer.flush();
					writer.close();
					
					randomize(unknownCausal, size, deg);
					randomize(unknownTarget, size, deg);
					
				} catch (IOException e) {
					e.printStackTrace();
				}
				count++;
				System.out.println(count+"/12 graphs generated.");
			}
		}
		System.out.println("Done making files!");
	}
	
	private class BloomFilter {
		private boolean[] bitList;
		
		public BloomFilter(int n) {
			bitList = new boolean[n*4];
			for(int i=0; i<bitList.length; i++) {
				bitList[i] = false;
			}
		}
		
		public void add(Edge e) {
			bitList[e.hash1() % bitList.length] = true;
			bitList[e.hash2() % bitList.length] = true;
			bitList[e.hash3() % bitList.length] = true;
			bitList[e.hash4() % bitList.length] = true;
			bitList[e.hash5() % bitList.length] = true;
			bitList[e.hash6() % bitList.length] = true;
			bitList[e.hash7() % bitList.length] = true;
			bitList[e.hash8() % bitList.length] = true;
//			bitList[e.hash9() % bitList.length] = true;
//			bitList[e.hash10() % bitList.length] = true;
		}
		
		public boolean contains(Edge e) {
			return bitList[e.hash1() % bitList.length]
					&& bitList[e.hash2() % bitList.length]
					&& bitList[e.hash3() % bitList.length]
					&& bitList[e.hash4() % bitList.length]
					&& bitList[e.hash5() % bitList.length]
					&& bitList[e.hash6() % bitList.length]
					&& bitList[e.hash7() % bitList.length]
					&& bitList[e.hash8() % bitList.length]
					/*&& bitList[e.hash9() % bitList.length]
					&& bitList[e.hash10() % bitList.length]*/;
		}
	}
	
	private class Edge {
		int head;
		int tail;
		
		public Edge(int h, int t) {
			head = h;
			tail = t;
		}
		
		@Override
		public boolean equals(Object o) {
			if(!(o instanceof Edge) || o == null) {
				return false;
			}
			Edge e = (Edge) o;
			return head == e.head && tail==e.tail;
		}
		
		@Override
		public int hashCode() {
			return head*tail;
		}
		
		public int hash1() {
			return head;
		}
		public int hash2() {
			return tail;
		}
		public int hash3() {
			return head + tail;
		}
		public int hash4() {
			return 2*head + tail;
		}
		public int hash5() {
			return head + 2*tail;
		}
		public int hash6() {
			return 3*head + 5*tail;
		}
		public int hash7() {
			return 7*head + 2*tail;
		}
		public int hash8() {
			return Math.abs((head-tail)*(head+tail));
		}
		public int hash9() {
			return 4*head + 11*tail;
		}
		public int hash10() {
			return Math.abs(5*head - 3*tail);
		}
	}
}
