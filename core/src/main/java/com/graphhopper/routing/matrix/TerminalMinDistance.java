package com.graphhopper.routing.matrix;

public class TerminalMinDistance {
    public int idx;
    public double weight;
    public long time;
    public double distance;

   public TerminalMinDistance(int idx,double weight, long time, double distance){
       this.idx = idx;
       this.weight = weight;
       this.time = time;
       this.distance = distance;
   }

    @Override
    public String toString() {
        return "TerminalMinDistance{" +
                "idx=" + idx +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                '}';
    }
}
