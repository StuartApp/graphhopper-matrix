package com.graphhopper.routing.matrix;

public class TurnCost {

    double weight;
    double cost;

    public TurnCost(double weight, double turncost){

        this.weight = weight;
        this.cost = turncost;
    }


    public double weightWithTurnCost(){
        return weight + cost;
    }

    public double turnCost(){
        return cost;
    }
}
