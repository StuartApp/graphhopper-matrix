package com.graphhopper.speeds;

public class Metrics {

    private double totalNumberEdgesProcessed;
    private double totalNumberEdgesConfiguredWithCustomSpeed;
    private double totalNumberEdgesConfiguredWithCustomSpeedInReverse;

    public Metrics() {
    }

    public double getTotalNumberEdgesProcessed() {
        return totalNumberEdgesProcessed;
    }

    public void incrementTotalNumberEdgesProcessed() {
        totalNumberEdgesProcessed = totalNumberEdgesProcessed + 1;
    }

    public double getTotalNumberEdgesConfiguredWithCustomSpeed() {
        return totalNumberEdgesConfiguredWithCustomSpeed;
    }

    public void incrementTotalNumberEdgesConfiguredWithCustomSpeed() {
        totalNumberEdgesConfiguredWithCustomSpeed = totalNumberEdgesConfiguredWithCustomSpeed + 1;
    }

    public double getTotalNumberEdgesConfiguredWithoutCustomSpeed() {
        return totalNumberEdgesProcessed - totalNumberEdgesConfiguredWithCustomSpeed;
    }

    public void incrementTotalNumberEdgesConfiguredWithCustomSpeedInReverse() {
        totalNumberEdgesConfiguredWithCustomSpeedInReverse = totalNumberEdgesConfiguredWithCustomSpeedInReverse + 1;
    }

    public double getTotalNumberEdgesConfiguredWithoutCustomSpeedInReverse() {
        return totalNumberEdgesProcessed - totalNumberEdgesConfiguredWithCustomSpeedInReverse;
    }
}
