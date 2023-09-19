/**
 * 
 * @author Nick FIsk
 *
 */
public class IndexPair {
	// a simple container class that holds 
	// the location of a seq in both the query and database arraylists.

	public int dbIndex;
	public int queIndex;
	
	public IndexPair(int dbIndex, int queIndex){
		this.dbIndex=dbIndex;
		this.queIndex=queIndex;
	}
	
	public int getDB(){
		return this.dbIndex;
	}
	public int getQUE(){
		return this.queIndex;
	}
}
