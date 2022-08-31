package com.graphhopper.routing.matrix.algorithm;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class SubGraph {

    protected List<Vertex> list[];

    public class Vertex {
        public int vertex;
        public double weight;
        public long time;
        public double distance;


        public Vertex(int vertex, double weight, long time, double distance) {
            this.vertex = vertex;
            this.weight = weight;
            this.time = time;
            this.distance = distance;
        }
    }

    int vertexes;

    public SubGraph(int vertexs) {
        this.vertexes = vertexs;
        list = new ArrayList[vertexs];
        for (int i = 0; i <vertexs ; i++) {
            list[i] = new ArrayList<>();
        }
    }

    public void addEdge(int source, int destination, double weight, long time, double distance){
        list[source].add(new Vertex(destination,weight,time,distance));
    }

    public List<Vertex> getNeighbours(int vertex){
        return list[vertex];
    }

    public void printGraph(){
        for (int i = 0; i < vertexes; i++) {
            if(list[i].size()>0) {
                System.out.print("Vertex " + i + " is connected to: ");
                for (int j = 0; j < list[i].size(); j++) {
                    System.out.print(list[i].get(j) + " ");
                }
                System.out.println();
            }
        }
    }
}
