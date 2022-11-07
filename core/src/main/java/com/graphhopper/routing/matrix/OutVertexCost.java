package com.graphhopper.routing.matrix;

public class OutVertexCost {

    TurnCost turnCost;
    Vertex out;


    public OutVertexCost(Vertex out, TurnCost turncost){

        this.out = out;
        this.turnCost = turncost;
    }


    public TurnCost getTurnCost() {
        return turnCost;
    }

    public Vertex getOut() {
        return out;
    }

    public double cost(){
        return this.turnCost.weightWithTurnCost() + out.weight;
    }
}
