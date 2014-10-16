package simulatedcount;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

public class SimulatedCount {

    public static void main(String[] args) {
        String dir="/home/jihongkim/workspace/chaoling/141016_simulation";
        String[] counts={"CL-10.bedtools.txt","CL-16.bedtools.txt","CL-22.bedtools.txt","CL-4.bedtools.txt","CL-12.bedtools.txt","CL-18.bedtools.txt"};
//        String[] counts={"CL-24.bedtools.txt","CL-6.bedtools.txt","CL-14.bedtools.txt","CL-20.bedtools.txt","CL-2.bedtools.txt","CL-7.bedtools.txt"};
        try {
            Map genes = new HashMap();
            FileReader fr = new FileReader(dir+"/simulated_annotation.txt");
            BufferedReader br = new BufferedReader(fr);
            String line="";
            while ((line=br.readLine())!=null) {
                String cols[]=line.split("\t");
                genes.put(cols[0], cols[1]);
            }
            br.close();
            fr.close();
            for (int i=0; i<counts.length; i++) {
                fr = new FileReader(dir+"/"+counts[i]);
                br = new BufferedReader(fr);
                FileWriter fw = new FileWriter(dir+"/"+counts[i].replaceAll(".txt", ".modified.txt"));
                BufferedWriter bw = new BufferedWriter(fw);
                while ((line=br.readLine())!=null) {
                    String cols[]=line.split("\t");
                    if (genes.get(cols[3])==null) bw.write(line+"\n");
                    else {
                        int cnt = Integer.parseInt(cols[12]);
                        if (genes.get(cols[3]).toString().equals("0")) cnt=(int) (cnt*((Math.random()*0.3d)+1d));
                    }
                }
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
