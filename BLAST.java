import java.io.BufferedReader;
//import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * 
 * @author Nick Fisk
 * BLAST is a rudimentary implementation of an older BLAST implementation. 
 * It finds hits to a query sequence for given data-sequences. 
 * Output is sent to standard out
 * 
 * Usage BLAST [FASTA Database] [FASTA query]
 *
 */
public class BLAST {
	//////////Timing things////////
	public static long overlapGen=0;
	public static long readData=0;
	public static long readQuery=0;
	public static long constructDatabase=0;
	public static long generateHits=0;
	public static long extendNtrim=0;
	public static long findUnique=0;
	public static long printData=0;
	public static long extend=0;
	public static long trim=0;
	public static long totalTime=0;
///////////////////////////////////////////////	
	
	
	public static final int wSize=28; //BLASTS standard sliding window size
	public static final int minSeed=5;//require at least 5 consequitive bases to start
	public static final String usage="usage: BLAST [FASTA database] [FASTA query]";//usage statement
	public static ArrayList<String> dbFASTAHeaders=new ArrayList<String>();//all the name info and such for the database FASTA
	public static ArrayList<String> dbSeqs=new ArrayList<String>(); //the actual sequences from database
	public static ArrayList<String> queryFASTAHeaders=new ArrayList<String>();//name info for query
	public static ArrayList<String> querySeqs=new ArrayList<String>();//sequence info for query
	////A hashmap representation of our database. The actual data is stored in the ArrayLists above
	////but the hashmap allows for quick lookup
	public static HashMap<String, ArrayList<IndexPair>>db=new HashMap<String, ArrayList<IndexPair>>();
	////hold all the hits
	public static ArrayList<Hit>hits=new ArrayList<Hit>();
	///hold only the unique hits, acts more as a checklist than anything. 
	public static ArrayList<Hit>singleHits=new ArrayList<Hit>();
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		////timin stuffs
		long totalTime1=System.nanoTime();
		long readData1=System.nanoTime();
		if(args.length!=2){//2 args
			System.out.println(usage);
			System.exit(1);
		}
		String databaseFile=args[0];
		String queryFile=args[1];
		ArrayList<String>allSeq=new ArrayList<String>();
		ArrayList<String>allSeq2=new ArrayList<String>();

		/////Read in the database
		
		BufferedReader read = new BufferedReader(
				new FileReader(databaseFile));
		String line = null;
		String curSeq="";
		
		
		while ((line = read.readLine()) != null) {
			allSeq.add(line.trim());
		}
		read.close();
		String holdMe="";
		for(String element: allSeq){
			if(element==""){
				System.out.println("WRONG");
			}
			if(element.charAt(0)=='>'){
				if(holdMe!=""){
					dbSeqs.add(holdMe);
					holdMe="";
				}
				dbFASTAHeaders.add(element.substring(1)); 
			}
			else{
				holdMe+=element;
			}
		}
		dbSeqs.add(holdMe);
		long readData2=System.nanoTime();
		readData+=readData2-readData1;
		long readQuery1=System.nanoTime();
		/////Read in the query
		read = new BufferedReader(
				new FileReader(queryFile));
		String line2 = null;
		//String curSeq2="";
		
		while ((line2 = read.readLine()) != null) {
			allSeq2.add(line2.trim());
		}
		String holdMe2="";
		for(String element: allSeq2){
			if(element==""){
				System.out.println("WRONG");
			}
			if(element.charAt(0)=='>'){
				if(holdMe2!=""){
					querySeqs.add(holdMe2);
					holdMe2="";
				}
				queryFASTAHeaders.add(element.substring(1)); 
			}
			else{
				holdMe2+=element;
			}
		}
		read.close();
		querySeqs.add(holdMe);
		long readQuery2=System.nanoTime();
		readQuery=readQuery2-readQuery1;
		long cdb1=System.nanoTime();
		
		///populate the database
		populateDatabase();
		long cdb2=System.nanoTime();
		constructDatabase=cdb2-cdb1;
		///timing is broken into pieces within findHits
		///find all the hits that meet the predefined criteria
		findHits();
		long totalTime2=System.nanoTime();
		totalTime=totalTime2-totalTime1;
		
		//////Print out the timing information. 
		System.out.println("\n\n\n\nRuntime Benchmarks");
		System.out.println("Database Read In: "+(readData/1000000)+ " milliseconds");
		System.out.println("Query Read In: "+(readQuery/1000000)+ " milliseconds");
		System.out.println("Database Construction: "+(constructDatabase/1000000)+ " milliseconds");
		System.out.println("Finding Hits: "+(generateHits/1000000)+ " milliseconds");
		System.out.println("Doing Hits and Overlap: "+(extendNtrim/1000000)+ " milliseconds");
		System.out.println("\t extending hits: " +(extend/1000000)+" milliseconds");
		System.out.println("\t  trimming hits: " +(trim/1000000)+" milliseconds");
		System.out.println("\t  generating overlap: " +(overlapGen/1000000)+" milliseconds");
		System.out.println("Finding unique hits: "+(findUnique/1000000)+ " milliseconds");
		System.out.println("Printing Data: "+(printData/1000000)+ " milliseconds (now part of unique hits)");
		System.out.println("Total Runtime: "+(totalTime/1000000)+ " milliseconds");
		
	}
	/**
	 * populates the database Hashmap. If a seed is already there, update the
	 * indecies it maps to. 
	 */
	public static void populateDatabase(){
		String sequence="";
		String chunk="";
		for(int i=0; i<dbSeqs.size(); i++){
			sequence=dbSeqs.get(i);
			for(int j=0; j<sequence.length()-wSize+1; j++){
				chunk=sequence.substring(j, j+wSize);
				if(db.containsKey(chunk)){
					ArrayList<IndexPair> temp=db.get(chunk);
					temp.add(new IndexPair(i,j));
					
				}
				else{
					ArrayList<IndexPair>temp= new ArrayList<IndexPair>();
					temp.add(new IndexPair(i,j));
					db.put(chunk, temp);
				}
			}
		}
	}
	/**
	 * Goes through the database and generates hits for each pairing.
	 * It also saves the unique sequences and prints out relevant informations
	 */
	public static void findHits(){
		long generateHits1=System.nanoTime();
		String sequence="";
		String chunk="";
		for(int i=0; i<querySeqs.size(); i++){
			sequence=querySeqs.get(i);
			for(int j=0; j<sequence.length()-wSize+1; j++){
				chunk=sequence.substring(j, j+wSize);
				if(db.containsKey(chunk)){
					ArrayList<IndexPair>temp=db.get(chunk);
					for(IndexPair p: temp){
						//System.out.println(p.dbIndex);
						//System.out.println(p.queIndex);
						//query, qheader, startquery, db, dbhead, start db
						//db, dbhead, query, qheader, start db, start query
						//System.out.println("ssize " +querySeqs.size());
						//System.out.println("hsize " +queryFASTAHeaders.size());
						//System.out.println(j);
						Hit htemp=(new Hit(dbSeqs.get(p.dbIndex),dbFASTAHeaders.get(p.dbIndex),sequence, queryFASTAHeaders.get(i) ,p.queIndex,j));
						hits.add(htemp);
						//System.out.println(htemp);
						//hits.add(new Hit(dbSeqs.get(temp.get(q).dbIndex),dbFASTAHeaders.get(temp.get(q).dbIndex),querySeqs.get(temp.get(q).queIndex), queryFASTAHeaders.get(temp.get(q).queIndex),i,j));
					}
				}
			}
		}
		long generateHits2=System.nanoTime();
		generateHits=generateHits2-generateHits1;
		long eat1=System.nanoTime();
		for(Hit h: hits){
			h.extendHit();
		}
		long eat2=System.nanoTime();
		extendNtrim=eat2-eat1;
		long fu1=System.nanoTime();
		int counter=1;
		for(Hit h: hits){
			if(!singleHits.contains(h)){
				h.mkMatch();
				singleHits.add(h);
				//System.out.println(h);
				System.out.print("///////////////////////");
				for(int i=0; i<Integer.toString(counter).length();i++){
					System.out.print("/");
				}
				System.out.println("///////////////////////");
				System.out.println("//////////  Hit Number "+counter+"  //////////////");
				System.out.print("///////////////////////");
				for(int i=0; i<Integer.toString(counter).length();i++){
					System.out.print("/");
				}
				System.out.println("///////////////////////");
				System.out.println(h);
				counter++;
			}
		}
		long fu2=System.nanoTime();
		findUnique=fu2-fu1;
		long pd1=System.nanoTime();
		/*int counter=1;
		for(Hit h: singleHits){
			h.mkMatch();
			System.out.print("///////////////////////");
			for(int i=0; i<Integer.toString(counter).length();i++){
				System.out.print("/");
			}
			System.out.println("///////////////////////");
			System.out.println("//////////  Hit Number "+counter+"  //////////////");
			System.out.print("///////////////////////");
			for(int i=0; i<Integer.toString(counter).length();i++){
				System.out.print("/");
			}
			System.out.println("///////////////////////");
			System.out.println(h);
			counter++;
		}*/
		long pd2=System.nanoTime();
		printData=pd2-pd1;
		
	}
	
}
