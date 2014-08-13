/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package com.lomolith.lncrna.finder;

import com.lomolith.common.model.Transcript;
import com.lomolith.sequence.LetterPairSimilarity;
import java.io.BufferedReader;
import java.io.BufferedWriter;
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
    
    ////////// DOUBLE CHECK FOR LENGTH-FILTER!!!!!
    ////////// DOUBLE CHECK FOR LIST ORDER!!!!!
    
    public static void main(String args[]) {
        String dir=".";
        String inputGTF="gencode.v19.lincRNA.gtf";
        String refGTF="gencode.v19.annotation.gtf";
        String refFasta="hg19.fa";
        int offset=0, range=100;
        String inputSeq="train.set";
        String output="output.txt";
        
        LncRNAFinder lncRNA = new LncRNAFinder(dir);
        boolean paired=true;
        boolean GET_SEQ=false;
        boolean GET_SIM=true;
        boolean GET_FREQ=false;
        int length=300;
        
        try {
            for (int i=0; i<args.length; i++) {
                if (args[i].equals("--seq")) GET_SEQ=true;                         // Get sequence of GTF files
                if (args[i].equals("--similarity")) GET_SIM=true;                           // Build a similarity score set
                if (args[i].equals("--count")) GET_FREQ=true;                           // Build a train set
                
                if (args[i].equals("-i") && i<args.length-1) inputGTF = args[i+1];
                if (args[i].equals("-d") && i<args.length-1) dir = args[i+1];
                if (args[i].equals("-g") && i<args.length-1) refGTF = args[i+1];
                if (args[i].equals("-f") && i<args.length-1) refFasta = args[i+1];
                if (args[i].equals("-s") && i<args.length-1) inputSeq = args[i+1];
                if (args[i].equals("-o") && i<args.length-1) output = args[i+1];
                if (args[i].equals("--offset") && i<args.length-1) offset = Integer.parseInt(args[i+1]);
                if (args[i].equals("--range") && i<args.length-1) range = Integer.parseInt(args[i+1]);
            }
            
            
            if (GET_SEQ) generateTrainSet(dir, inputGTF, refFasta, lncRNA, inputSeq);
            else if (GET_FREQ) {
                /// CHECK TRAINSET FIRST!!
                List trainset = readTrainSequence(dir, inputSeq);
//                for (int i=offset*range; i<(offset+1)*range; i++) {
                String nuc[]={"A","C","T","G"};
                FileWriter fw = new FileWriter(dir+"/"+output);
                BufferedWriter bw = new BufferedWriter(fw);
                bw.write("ID\tLength\tORF length\tORF\tA\tT\tC\tG");
                for (int j1=0; j1<4; j1++)
                    for (int j2=0; j2<4; j2++)
                        for (int j3=0; j3<4; j3++) {
                            String idx=nuc[j1]+nuc[j2]+nuc[j3];
                            bw.write("\t"+idx);
                        }
                bw.write("\n");
                    
                for (int i=0; i<trainset.size(); i++) {
                    Transcript tc1=(Transcript) trainset.get(i);
                    System.out.println("Counting tri-nucleotides for "+tc1.id+" ("+(i+1)+"/"+trainset.size()+")...");
                    String seq=tc1.sequence.toUpperCase();
                    
                    boolean start=false;
                    if (seq.indexOf("ATG")!=-1) {start=true; seq=seq.substring(seq.indexOf("ATG")+3);}

//                    int ochre=seq.indexOf("TAA");
//                    int amber=seq.indexOf("TAG");
//                    int opal_umber=seq.indexOf("TGA");
//                    if (ochre!=-1 && amber==-1 && opal_umber==-1) seq=seq.substring(0,ochre);
//                    else if (ochre==-1 && amber!=-1 && opal_umber==-1) seq=seq.substring(0,amber);
//                    else if (ochre==-1 && amber==-1 && opal_umber!=-1) seq=seq.substring(0,opal_umber);
//                    else if (ochre==-1 && amber!=-1 && opal_umber!=-1) seq=seq.substring(0,Math.min(amber, opal_umber));
//                    else if (ochre!=-1 && amber==-1 && opal_umber!=-1) seq=seq.substring(0,Math.min(ochre, opal_umber));
//                    else if (ochre!=-1 && amber!=-1 && opal_umber==-1) seq=seq.substring(0,Math.min(ochre, amber));
//                    else if (ochre!=-1 && amber!=-1 && opal_umber!=-1) seq=seq.substring(0,Math.min(ochre, Math.min(amber, opal_umber)));
                    Map cnt = new HashMap();
                    boolean closed=false;
                    int pos=0;                    
                    for (int j=0; j<seq.length()-3; j+=3) {
                        String idx = seq.substring(j,j+3);
                        if (idx.equals("TAA")||idx.equals("TAG")||idx.equals("TGA")) {closed=true; pos=j; break;}
                        if (cnt.get(idx)==null) cnt.put(idx, "0");
                        int c = Integer.parseInt(cnt.get(idx).toString())+1;
                        cnt.put(idx, c);
                    }
//System.out.println(seq.substring(0,pos+3));
                    bw.write(tc1.id+"\t"+tc1.sequence.length()+"\t"+pos);
                    if (start && closed) bw.write("\tYES");
                    else if (start && !closed) bw.write("\tSTART");
                    else if (!start && closed) bw.write("\tSTOP");
                    else bw.write("\tNO");

                    String single_freq=seq;                                     // Count only in ORF
                    String remove_A=single_freq.replaceAll("A", "");
                    String remove_T=remove_A.replaceAll("T", "");
                    String remove_C=remove_T.replaceAll("C", "");
                    String remove_G=remove_C.replaceAll("G", "");
                    bw.write("\t"+(single_freq.length()-remove_A.length())+"\t"+(remove_A.length()-remove_T.length())+"\t"+(remove_T.length()-remove_C.length())+"\t"+(remove_C.length()-remove_G.length()));
                    
                    for (int j1=0; j1<4; j1++)
                        for (int j2=0; j2<4; j2++)
                            for (int j3=0; j3<4; j3++) {
                                String idx=nuc[j1]+nuc[j2]+nuc[j3];
                                bw.write("\t"+(cnt.get(idx)==null?"0":cnt.get(idx)));
                            }
                    bw.write("\n");
                }
                bw.close();
                fw.close();
            }
            else if (GET_SIM) {
                List trainset = readTrainSequence(dir, inputSeq);
                for (int i=offset*range; i<(offset+1)*range; i++) {
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
        List trainset = lncRNA.getTrainSequences(dir+"/"+inputGTF, refFasta, false);
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
//        List trainset = readTrainSequence(dir);
        
//        GTF ref = readRef(dir, "gencode.v19.annotation.gtf");
//        ref.sortTranscripts();
//        List ref_trans = ref.getTranscripts();
//        
//        FileWriter fw = new FileWriter(dir+"/train.set.freq_dist");
//        BufferedWriter bw = new BufferedWriter(fw);
//        bw.write("Length\tA\tT\tC\tG\tDist\n");
//        for (int i=0; i<trainset.size(); i++) {
//            Transcript tc = (Transcript) trainset.get(i);
//            String seq_A = tc.sequence.toUpperCase().replaceAll("A", "");
//            String seq_T = seq_A.replaceAll("T","");
//            String seq_C = seq_T.replaceAll("C","");
//            long proximity=9999999;
//            for (int j=0; j<ref_trans.size(); j++) {
//                Transcript t = (Transcript) ref_trans.get(j);
//                String chr1=tc.chr.toUpperCase().replaceAll("CHR","");
//                String chr2=t.chr.toUpperCase().replaceAll("CHR","");
//                if (chr1.equals("X") && !chr2.equals("Y") || chr1.equals("Y") || (!chr2.equals("X") && !chr2.equals("Y") && Integer.parseInt(chr1)<=Integer.parseInt(chr2))) {
//                    if (chr1.equals(chr2)) {
//                        if (tc.start<=t.end && tc.end>=t.start) proximity=-1;
//                        else proximity=Math.min(Math.abs(tc.start-t.end), Math.abs(tc.end-t.start));
//                    }
//                }
//                else break;
//            }
//            bw.write(tc.sequence.length()+"\t"+(tc.sequence.length()-seq_A.length())+"\t"+(seq_A.length()-seq_T.length())+"\t"+(seq_T.length()-seq_C.length())+"\t"+seq_C.replaceAll("T","").length()+"\t"+proximity+"\n");
//        }
//        bw.close();
//        fw.close();
//        
//        
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
