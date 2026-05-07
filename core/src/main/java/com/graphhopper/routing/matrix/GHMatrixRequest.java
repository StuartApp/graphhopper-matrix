package com.graphhopper.routing.matrix;

import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.graphhopper.util.CustomModel;
import com.graphhopper.util.Helper;
import com.graphhopper.util.PMap;
import com.graphhopper.util.Parameters;
import com.graphhopper.util.shapes.GHPoint;

import java.util.ArrayList;
import java.util.List;

public class GHMatrixRequest {
    // GHMatrixRequest doesn't have any hint because it is just possible to work just for now with CH, so the profile
    // settings cannot be changed. On the other hand, for now we don't want to open the API to be configurable, to reduce the complexity
    // of the flow.

    private String algo = Parameters.Algorithms.DIJKSTRA_MANY_TO_MANY;
    private String profile = "";
    private final List<GHPoint> origins = new ArrayList<>();
    private final List<GHPoint> destinations = new ArrayList<>();
    private final List<String> originPointHints = new ArrayList<>();
    private final List<String> destinationPointHints = new ArrayList<>();
    private List<String> snapPreventions = new ArrayList<>();
    private List<String> originCurbsides = new ArrayList<>();
    private List<String> destinationCurbsides = new ArrayList<>();
    private final PMap hints = new PMap();
    private boolean failFast = true;
    private CustomModel customModel;

    /**
     * One or more locations to use as the starting point for calculating travel distance and time.
     */
    public List<GHPoint> getOrigins() {
        return origins;
    }

    @JsonProperty("from_points")
    public GHMatrixRequest setOrigins(List<GHPoint> origins) {
        this.origins.addAll(origins);
        return this;
    }

    /**
     * One or more locations to use as the finishing point for calculating travel distance and time.
     */
    public List<GHPoint> getDestinations() {
        return destinations;
    }

    @JsonProperty("to_points")
    public GHMatrixRequest setDestinations(List<GHPoint> destinations) {
        this.destinations.addAll(destinations);
        return this;
    }

    public GHMatrixRequest setPoints(List<GHPoint> origins) {
        this.setOrigins(origins);
        this.setDestinations(origins);
        return this;
    }

    public String getAlgorithm() {
        return algo;
    }

    public GHMatrixRequest setAlgorithm(String algo) {
        if (algo != null)
            this.algo = Helper.camelCaseToUnderScore(algo);
        return this;
    }

    public String getProfile() {
        return profile;
    }

    public GHMatrixRequest setProfile(String profile) {
        this.profile = profile;
        return this;
    }

    public PMap getHints() {
        return hints;
    }

    /**
     * This method sets a key value pair in the hints and is unrelated to the setPointHints method.
     * It is mainly used for deserialization with Jackson.
     *
     * @see #setPointHints(List)
     */
    // a good trick to serialize unknown properties into the HintsMap
    @JsonAnySetter
    public GHMatrixRequest putHint(String fieldName, Object value) {
        this.hints.putObject(fieldName, value);
        return this;
    }

    public boolean isFailFast() {
        return failFast;
    }

    public GHMatrixRequest setFailFast(boolean failFast) {
        this.failFast = failFast;
        return this;
    }

    public List<String> getOriginPointHints() {
        return originPointHints;
    }

    @JsonProperty("from_point_hints")
    public GHMatrixRequest setOriginPointHints(List<String> pointHints) {
        this.originPointHints.addAll(pointHints);
        return this;
    }

    public List<String> getDestinationPointHints() {
        return destinationPointHints;
    }

    @JsonProperty("to_point_hints")
    public GHMatrixRequest setDestinationPointHints(List<String> pointHints) {
        this.destinationPointHints.addAll(pointHints);
        return this;
    }

    @JsonProperty("point_hints")
    public GHMatrixRequest setPointHints(List<String> pointHints) {
        this.originPointHints.addAll(pointHints);
        this.destinationPointHints.addAll(pointHints);
        return this;
    }

    public List<String> getSnapPreventions() {
        return snapPreventions;
    }

    public GHMatrixRequest setSnapPreventions(List<String> snapPreventions) {
        this.snapPreventions.addAll(snapPreventions);
        return this;
    }

    public List<String> getOriginCurbsides() {
        return originCurbsides;
    }

    @JsonProperty("from_curbsides")
    public GHMatrixRequest setOriginCurbsides(List<String> curbsides) {
        this.originCurbsides.addAll(curbsides);
        return this;
    }

    public List<String> getDestinationCurbsides() {
        return destinationCurbsides;
    }

    @JsonProperty("to_curbside")
    public GHMatrixRequest setDestinationCurbsides(List<String> curbsides) {
        this.destinationCurbsides.addAll(curbsides);
        return this;
    }

    public GHMatrixRequest setCurbsides(List<String> curbsides) {
        this.destinationCurbsides.addAll(curbsides);
        this.originCurbsides.addAll(curbsides);
        return this;
    }

    public CustomModel getCustomModel() {
        return customModel;
    }

    public GHMatrixRequest setCustomModel(CustomModel customModel) {
        this.customModel = customModel;
        return this;
    }
}
