package com.graphhopper.routing.matrix;

import java.util.Objects;

public class Terminal {

    public double weight;
    public long time;
    public double distance;
    public int node;
    public int nodeIdx;

    public Terminal(Vertex vert, int terminal, int terminalIdx) {
        this.weight = vert.weight;
        this.time = vert.time;
        this.distance = vert.distance;
        this.node = terminal;
        this.nodeIdx = terminalIdx;

    }

    protected Terminal(double weight, long time, double distance, int terminal, int terminalIdx) {
        this.weight =weight;
        this.time = time;
        this.distance = distance;
        this.node = terminal;
        this.nodeIdx = terminalIdx;

    }

    public Terminal with(double weight, long time, double distance){
        return new Terminal(this.weight + weight, this.time + time, this.distance + distance,
                this.node,this.nodeIdx);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Terminal terminal = (Terminal) o;
        return nodeIdx == terminal.nodeIdx;
    }

    @Override
    public int hashCode() {
        return Objects.hash(nodeIdx);
    }


    @Override
    public String toString() {
        return "Terminal{" +
                "weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                ", terminal=" + node +
                ", terminalIdx=" + nodeIdx +
                '}';
    }
}
