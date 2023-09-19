/*
 * Hit represents a hit to the database
 * The toString is actually how the print in BLAST is formatted.
 * It has a slew of accessor functions, few of which are used in this version.
 */
/**
 * 
 * @author Nick Fisk
 *
 */
public class Hit {
	public String dbSeq;
	public String dbFASTA;
	public String querySeq;
	public String queryFASTA;
	public int qStart;
	public int qEnd;
	public int dStart;
	public int dEnd;
	public String Match;//the overlap
	public static final int allowableMismatches=5;//the number of mismatches allowed as we extend.
	public static int count=0;
	public static int count2=0;
	//public static int timesCalled=0;
	
	/*
	 * Constructor for a Hit. Pretty much toss all relevant information about the hit into the object for easy
	 * refference.
	 */
	public Hit(String dbSeq,String dbFASTA, String querySeq, String queryFASTA, int dStart, int qStart){
		this.dbSeq=dbSeq;
		this.dbFASTA=dbFASTA;
		this.querySeq=querySeq;
		this.queryFASTA=queryFASTA;
		this.qStart=qStart;
		this.dStart=dStart;
		this.qEnd=qStart+BLAST.wSize;
		this.dEnd=dStart+BLAST.wSize;
		//if(count<5){
		//System.out.println(qStart+" "+ dStart);
		//System.out.println(dbSeq.length()+" "+querySeq.length());
		//count++;}
		
		
	}
	
	/////Accessor Functions out da wazoo, should be an easy optimization////
	
	public String getdbSeq(){
		return this.dbSeq;
	}
	public String getdbFASTA(){
		return this.dbFASTA;
	}
	public String getquerySeq(){
		return this.querySeq;
	}
	public String getqueryFASTA(){
		return this.queryFASTA;
	}
	public int getqStart(){
		return this.qStart;
	}
	public int getqEnd(){
		return this.qEnd;
	}
	public int getdStart(){
		return this.dStart;
	}
	public int getdEnd(){
		return this.dEnd;
	}
	public String getMatch(){
		return this.Match;
	}
	
	/////Printable format of a hit
	@Override
	public String toString(){
		String toPrint="";
		int index;
		int maxLen=dEnd-dStart;
		System.out.print(qStart+" "+qEnd+"\n");
		for(int i=0; i<maxLen; i+=57){
			if(i+57>maxLen){
				index=maxLen; ///just like printing the MSA
			}
			else{
				index=i+57;
			}
			///Amanda suggests using String.format. Give it a go!
			toPrint+=String.format("%-10s%5s ", dbFASTA, dStart+i);
			toPrint+=(dbSeq.substring(dStart+i,dStart+index));
			toPrint+=(String.format(" %5s", dStart+index));
			toPrint+="\n";//start query info
			toPrint+=String.format("%-10s%5s ", queryFASTA, qStart+i);
			toPrint+=querySeq.substring(qStart+i, qStart+index);
			toPrint+=String.format(" %5s", qStart+index);
			toPrint+="\n"; //start overlap info
			toPrint+=String.format("%-16s", "");
			toPrint+=Match.substring(i, index);
			toPrint+="\n\n\n";
				
		}
		
		return toPrint;
		
	}
	///can't use getOverlap----is the name of the accesor.
	/**
	 * mkMatch makes the overlap sequence (kind of like a consensus sequence, but simpler.
	 */
	public void mkMatch(){
		long t1=System.nanoTime();
		String overlap="";
		for(int i=0; i<getdEnd()-dStart; i++){
			if(querySeq.charAt(qStart+i)==dbSeq.charAt(dStart+i)){
				overlap+="*";
			}
			else{
				overlap+="_";
			}
		}
		this.Match=overlap;
		long t2=System.nanoTime();
		BLAST.overlapGen+=(t2-t1);
	}
	//extend the left side of the seed until a mismatch threshhold is exceeded
	public void doLeft(){

		int leftMM=0;
		
		//timesCalled+=1;
		//System.out.println(timesCalled);
		while(leftMM<allowableMismatches&&0<qStart&&0<dStart){
			if((querySeq.charAt(qStart-1))==dbSeq.charAt(dStart-1)){
				leftMM=0;
			}
			else{
				leftMM+=1;
			}
			dStart--;
			qStart--;
		}
	}
	
	////extend the right seed until mismatch threshold
	public void doRight(){
		
		int rightMM=0;
		while(rightMM<allowableMismatches&&dEnd<dbSeq.length()&&qEnd<querySeq.length()){
			if(querySeq.charAt(qEnd)==dbSeq.charAt(dEnd)){
				rightMM=0;
			}
			else{
				rightMM+=1;
			}
			this.dEnd++;
			this.qEnd++;
		}
		
	}
	
	///trim back the seed (index) to the appropriate location
	public void trim(){
	//	if(count2<5){
//			count2++;
		//System.out.println(dStart+", " +qStart+ ", " +dEnd+", "+ qEnd);
	//}
		long t1=System.nanoTime();
		while(querySeq.charAt(qStart)!=dbSeq.charAt(dStart)){
			this.dStart+=1;
			this.qStart+=1;
			System.out.println(this.qStart + "  ");
		}
		while(querySeq.charAt(qEnd-1)!=dbSeq.charAt(dEnd-1)){
			this.dEnd-=1;
			this.qEnd-=1;
			if(qEnd-1==-1||dEnd-1==-1){
				System.out.println("why");
				break;
			}
		}
		long t2=System.nanoTime();
		BLAST.trim+=(t2-t1);
	}
	///called by BLAST to perform the major functions of Hit
	public void extendHit(){
		long t1=System.nanoTime();
		doLeft();
		doRight();
		long t2=System.nanoTime();
		BLAST.extend+=t2-t1;
		//doLeft();
		trim();
		//mkMatch();
	}
	/////This is needed to use the .contains method of the arrayList XD 
	@Override
	public boolean equals(Object o){
		boolean areEqual=false;
		if(o!=null){
			if(o instanceof Hit){
				Hit h=(Hit)o;
				areEqual=(this.querySeq==h.getquerySeq()&&
						this.dbSeq==h.getdbSeq()&&
						this.dStart==h.getdStart()&&
						this.qStart==h.getqStart()&&
						this.dEnd==h.getdEnd()&&
						this.qEnd==h.getqEnd()
						);
				/*if((this.querySeq.equals(h.getquerySeq()))&&this.dbSeq.equals(h.getdbSeq())&&
						(this.dEnd==h.getdEnd())&&(this.qEnd==h.getqEnd())&&(this.dStart
						==h.getdStart())&&(this.qStart==h.getqStart())&&(this.dStart==getqEnd())){
					areEqual=true;
				}*/
			}
		}
		
		return areEqual;
		
	}
}
