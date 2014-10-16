package simulatedcount;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class SimulatedCount {

    public static void main(String[] args) {
        String dir="/home/jihongkim/workspace/chaoling/141016_simulation";
        String[] counts={"CL-10.bedtools.txt","CL-16.bedtools.txt","CL-22.bedtools.txt","CL-4.bedtools.txt","CL-12.bedtools.txt","CL-18.bedtools.txt","CL-24.bedtools.txt","CL-6.bedtools.txt","CL-14.bedtools.txt","CL-20.bedtools.txt","CL-2.bedtools.txt","CL-7.bedtools.txt"};
        try {
            for (int i=0; i<counts.length; i++) {
                FileReader fr = new FileReader(dir+"/"+counts[i]);
                BufferedReader br = new BufferedReader(fr);
                FileWriter fw = new FileWriter(dir+"/"+counts[i].replaceAll(".txt", ".modified.txt"));
                BufferedWriter bw = new BufferedWriter(fw);
                
                
                bw.close();
                fw.close();
                br.close();
                fr.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
