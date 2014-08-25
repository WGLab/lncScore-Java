package com.lomolith.lncrna.finder;

import com.lomolith.common.model.Transcript;
import com.lomolith.sequence.LetterPairSimilarity;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author mittens
 */
public class Runner {
    
    static boolean VERBOSE=true;
    
    // REMOVE DUPLICATES IN DATA
    // REMOVE OVERLAP IN TEST DATA
    
    // Possible GUI filter option: Length
    //                             ORF closed? exist? length?
    
    
    public static void main(String args[]) {
        String dir=".";
        String inputGTF="gencode.v19.lincRNA.gtf";
        String refGTF="gencode.v19.annotation.gtf";
        String refFasta="hg19.fa";
        int subset=0, range=100, offset=6;
        String inputSeq="train.set";
        String output="output.txt";
        String head_file="";

        String fileCDS="";
        String fileMEF="";
        String filePhastCon="";
        String filelincRNA="";
        
        LncRNAFinder lncRNA = new LncRNAFinder(dir);
        boolean GET_SEQ=false;
        boolean GET_SIM=false;
        boolean GET_COUNT=false;
        boolean GET_FREQ=false;
        boolean GET_INPUT=false;
        boolean FORMAT_LIBSVM=false;
        int length=300;
        try {
            for (int i=0; i<args.length; i++) {
                if (args[i].equals("--seq")) GET_SEQ=true;                         // Get sequence of GTF files
                if (args[i].equals("--similarity")) GET_SIM=true;                           // Build a similarity score set
                if (args[i].equals("--count")) GET_COUNT=true;                           // Build a train set
                if (args[i].equals("--freq")) GET_FREQ=true;
                if (args[i].equals("--generate")) GET_INPUT=true;                           // Build a similarity score set
                
                if (args[i].equals("-i") && i<args.length-1) inputGTF = args[i+1];
                if (args[i].equals("-d") && i<args.length-1) dir = args[i+1];
                if (args[i].equals("-g") && i<args.length-1) refGTF = args[i+1];
                if (args[i].equals("-f") && i<args.length-1) refFasta = args[i+1];
                if (args[i].equals("-s") && i<args.length-1) inputSeq = args[i+1];
                if (args[i].equals("-o") && i<args.length-1) output = args[i+1];
                if (args[i].equals("--subset") && i<args.length-1) subset = Integer.parseInt(args[i+1]);
                if (args[i].equals("--libsvm")) FORMAT_LIBSVM=true;
                if (args[i].equals("--offset") && i<args.length-1) offset = Integer.parseInt(args[i+1]);
                if (args[i].equals("--header") && i<args.length-1) head_file = args[i+1];
                if (args[i].equals("--range") && i<args.length-1) range = Integer.parseInt(args[i+1]);
                
                if (args[i].equals("--cds") && i<args.length-1) fileCDS = args[i+1];
                if (args[i].equals("--phastcon") && i<args.length-1) filePhastCon = args[i+1];
                if (args[i].equals("--mef") && i<args.length-1) fileMEF = args[i+1];
                if (args[i].equals("--lincRNA") && i<args.length-1) filelincRNA = args[i+1];
            }
            
            if (GET_SEQ) generateTrainSet(dir, inputGTF, refFasta, lncRNA, inputSeq);
            else if (GET_COUNT || GET_FREQ || GET_INPUT) {
                /// CHECK TRAINSET FIRST!!
                Map feat = new HashMap();
                if (GET_INPUT) {
                    System.out.println("Reading lincRNA list...");
                    FileReader fr = new FileReader(dir+"/"+filelincRNA);
                    BufferedReader br = new BufferedReader(fr);
                    String line="";
                    while((line=br.readLine())!=null) {
                        String cols[]=line.split("\t");
                        feat.put(cols[0],cols[1]);
                    }
                    System.out.println("Importing coding potential...");
                    fr = new FileReader(dir+"/"+fileCDS);
                    br = new BufferedReader(fr);
                    while((line=br.readLine())!=null) {
                        String cols[]=line.split("\t");
                        if (feat.get(cols[0])!=null) {
                            String out=feat.get(cols[0]).toString()+"\t1:"+cols[5];
                            feat.put(cols[0],out);
                        }
                    }
                    System.out.println("Importing conservation score...");
                    fr = new FileReader(dir+"/"+filePhastCon);
                    br = new BufferedReader(fr);
                    while((line=br.readLine())!=null) {
                        String cols[]=line.split("\t");
                        if (feat.get(cols[7])!=null) {
                            String out=feat.get(cols[7]).toString()+"\t2:"+cols[1].substring(cols[1].indexOf("=")+1,cols[1].indexOf(";"))+"\t3:"+cols[1].substring(cols[1].lastIndexOf("=")+1);
                            feat.put(cols[7],out);
                        }
                    }
                    System.out.println("Importing minimal folding energy...");
                    fr = new FileReader(dir+"/"+fileMEF);
                    br = new BufferedReader(fr);
                    while((line=br.readLine())!=null) {
                        String id=line.trim().substring(1);
                        line=br.readLine();
                        if (line.startsWith("--")) line=br.readLine();
                        line=line.trim().replaceAll("[()]", "");
                        if (feat.get(id)!=null) {
                            String out=feat.get(id).toString()+"\t4:"+line;
                            feat.put(id,out);
System.out.println(id+"\t"+out);                        
                        }
                    }
                    br.close();
                    fr.close();
                    offset=5;
                }
                else if (!head_file.equals("") && (new File(dir+"/"+head_file)).exists()) {
                    System.out.println("Importing exist features...");
                    FileReader fr = new FileReader(dir+"/"+head_file);
                    BufferedReader br = new BufferedReader(fr);
                    String line=br.readLine();
                    String cols[]=line.split("\t");
                    offset=cols.length-1;
                    while((line=br.readLine())!=null) {
                        cols=line.split("\t");
                        String out=cols[1];
                        for (int i=2; i<cols.length; i++) out+="\t"+(i-1)+":"+cols[i];
                        feat.put(cols[0],out);
                    }
                    br.close();
                    fr.close();
                }
                
                List trainset = readTrainSequence(dir, inputSeq);
                String nuc[]={"A","C","T","G"};
                FileWriter fw = new FileWriter(dir+"/"+output);
                BufferedWriter bw = new BufferedWriter(fw);
                if (!FORMAT_LIBSVM) {
                    bw.write("ID\tLength\tORF length\tORF\tA\tT\tC\tG");
                    for (int j1=0; j1<4; j1++)
                        for (int j2=0; j2<4; j2++)
                            for (int j3=0; j3<4; j3++) {
                                String idx=nuc[j1]+nuc[j2]+nuc[j3];
                                bw.write("\t"+idx);
                            }
                    bw.write("\n");
                }
                    
                for (int i=0; i<trainset.size(); i++) {
                    Transcript tc1=(Transcript) trainset.get(i);
                    if (head_file.equals("") || feat.get(tc1.id)!=null) {
                        System.out.println("Counting tri-nucleotides for "+tc1.id+" ("+(i+1)+"/"+trainset.size()+")...");
                        String seq=tc1.sequence.toUpperCase();

                        boolean start=false;
                        if (seq.indexOf("ATG")!=-1) {start=true; seq=seq.substring(seq.indexOf("ATG")+3);}

                        Map cnt = new HashMap();
                        boolean closed=false;
                        int pos=0;                    
                        int total_cnt=0;
                        for (int j=0; j<seq.length()-3; j+=3) {
                            String idx = seq.substring(j,j+3);
                            if (idx.equals("TAA")||idx.equals("TAG")||idx.equals("TGA")) {closed=true; pos=j;}
                            else {if (cnt.get(idx)==null) cnt.put(idx, "0");
                                int c = Integer.parseInt(cnt.get(idx).toString())+1;
                                cnt.put(idx, c);
                                total_cnt++;
                            }
                        }
                        int index=offset;
                        if (!FORMAT_LIBSVM) {
                            bw.write(tc1.id+"\t"+tc1.sequence.length()+"\t"+pos);
                            if (start && closed) bw.write("\tYES");
                            else if (start && !closed) bw.write("\tSTART");
                            else if (!start && closed) bw.write("\tSTOP");
                            else bw.write("\tNO");
                            bw.write(tc1.id+"\t"+(index++)+":"+tc1.sequence.length()+"\t"+(index++)+":"+pos); 
                        }
                        else {
                            if (!head_file.equals("")||GET_INPUT) bw.write(feat.get(tc1.id).toString()+"\t");
                            bw.write((index++)+":"+tc1.sequence.length()+"\t"+(index++)+":"+pos); 
                        }

                        String single_freq=seq;                                     // Count only in ORF
                        String remove_A=single_freq.replaceAll("A", "");
                        String remove_T=remove_A.replaceAll("T", "");
                        String remove_C=remove_T.replaceAll("C", "");
                        String remove_G=remove_C.replaceAll("G", "");
                        if (GET_COUNT) bw.write("\t"+((single_freq.length()-remove_A.length())/seq.length())+"\t"+((remove_A.length()-remove_T.length())/seq.length())+"\t"+(remove_T.length()-remove_C.length())+"\t"+(remove_C.length()-remove_G.length()));
                        else bw.write("\t"+(index++)+":"+((double)(single_freq.length()-remove_A.length())/(double)seq.length())+"\t"+(index++)+":"+((double)(remove_A.length()-remove_T.length())/(double)seq.length())+"\t"+(index++)+":"+((double)(remove_T.length()-remove_C.length())/(double)seq.length())+"\t"+(index++)+":"+((double)(remove_C.length()-remove_G.length())/(double)seq.length()));

                        for (int j1=0; j1<4; j1++)
                            for (int j2=0; j2<4; j2++)
                                for (int j3=0; j3<4; j3++) {
                                    String idx=nuc[j1]+nuc[j2]+nuc[j3];
                                    if (GET_COUNT) bw.write("\t"+(cnt.get(idx)==null?"0":cnt.get(idx)));
                                    else bw.write("\t"+(index++)+":"+(cnt.get(idx)==null?"0":(Double.parseDouble(cnt.get(idx).toString())/((double)total_cnt))));
                                }
                        bw.write("\n");
                    }
                }
                bw.close();
                fw.close();
            }
            else if (GET_SIM) {
                List trainset = readTrainSequence(dir, inputSeq);
                for (int i=subset*range; i<(subset+1)*range; i++) {
                    System.out.println("Calculating simirarity ("+(i+1)+"/"+trainset.size()+")...");
                    Transcript tc1=(Transcript) trainset.get(i);
                    double len=0.0;
                    double min=9999;
                    double max=0;
                    double similarity=0;
                    for (int j=0; j<trainset.size(); j++) {
                        if (i!=j) {
                            Transcript tc2=(Transcript) trainset.get(j);
                            len=len+tc2.sequence.length();
                            double sim=LetterPairSimilarity.compareStrings(tc1.sequence, tc2.sequence);
                            similarity+=sim;
                            if (min>sim) min=sim;
                            if (max<sim) max=sim;
                        }
                    }
                    System.out.println("   Length:"+tc1.sequence.length()+"("+(len/(trainset.size()-1))+")\tSimilarity:"+(similarity/(trainset.size()-1))+"("+min+"~"+max+")");
                }
            }
            else {  /// No option
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void generateTrainSet(String dir, String inputGTF, String refFasta, LncRNAFinder lncRNA, String inputSeq) throws IOException {
        lncRNA = new LncRNAFinder(dir);
        List trainset = lncRNA.getTrainSequences(inputGTF, refFasta, false);
        FileWriter fw = new FileWriter(dir+"/"+inputSeq);
        FileWriter fw2 = new FileWriter(dir+"/"+inputSeq+".fa");
        BufferedWriter bw = new BufferedWriter(fw);
        BufferedWriter bw2 = new BufferedWriter(fw2);
        for (Object t: trainset) {
            Transcript tc=(Transcript) t;
            bw.write(tc.id+"\t"+tc.sequence+"\n");
            bw2.write(">"+tc.id+"\n");
            String output=tc.sequence;
            while(output.length()>80) {
                bw2.write(output.substring(0,80)+"\n");
                output=output.substring(80);
            }
            bw2.write(output+"\n");
        }
        bw2.close();
        fw2.close();
        bw.close();
        fw.close();
    }
    
    private static List readTrainSequence(String dir, String seqeunceFile) throws FileNotFoundException, IOException {
        List trainset = new ArrayList();
        FileReader fr = new FileReader(dir+"/"+seqeunceFile);
        BufferedReader br = new BufferedReader(fr);
        String line="";
        while ((line=br.readLine())!=null) {
            Transcript t = new Transcript();
            String cols[]=line.split("\t");
            t.id=cols[0];
            t.sequence=cols[1];
            trainset.add(t);
        }
        br.close();
        fr.close();
        return trainset;
    }
}
