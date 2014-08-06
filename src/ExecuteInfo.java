
public class ExecuteInfo {
	double time;
	int numItr;
	long maxHeap, endHeap;
	
	/**
	 * @purpose Used for flowStart() in SingleSourceKShortestPaths.
	 */
	ExecuteInfo(double t, long eHeap, long mHeap) {
		time = t;
		endHeap = eHeap;
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
}
