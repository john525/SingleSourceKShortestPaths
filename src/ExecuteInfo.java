
public class ExecuteInfo {
	double time;
	int numItr;
	long maxHeap;
	//double timeToReadMem;
	
	/**
	 * @purpose Used for flowStart() in SingleSourceKShortestPaths.
	 */
	ExecuteInfo(double t, long mHeap) {
		time = t;
		maxHeap = mHeap;
	}
	
	/**
	 * @purpose Used for flowStart() in SingleSourceKDiverseShortPaths.
	 */
	ExecuteInfo(double t) {
		time = t;
	}
	
	/**
	 * @purpose Used by SingleSourceKDiverseShortPaths multipleRuns().
	 */
	ExecuteInfo(double t, int itr){
		this.time = t;
		numItr = itr;
	}
	
	/**
	 * 
	 * @return The number of bytes in the heap. Should be normalized by subtracting
	 * its original value from the start of the program.
	 */
	public static long memoryUsed() {
		Runtime rt = Runtime.getRuntime();
		rt.gc();
		return rt.totalMemory() - rt.freeMemory();
	}
	
	public String toString() {
		return "Time: " + time  +"s\nMemory:"+(maxHeap/1e6)+"MB";
	}
}
