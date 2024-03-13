package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.*;
import com.graphhopper.storage.IntsRef;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class OSMTollMotorcycleParserTest {
    private EnumEncodedValue<TollMotorcycle> tollEnc;
    private OSMTollMotorcycleParser parser;

    @BeforeEach
    public void setUp() {
        tollEnc = new EnumEncodedValue<>(TollMotorcycle.KEY, TollMotorcycle.class);
        tollEnc.init(new EncodedValue.InitializerConfig());
        parser = new OSMTollMotorcycleParser(tollEnc);
    }

    @Test
    public void testTollMotorcycleYes() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:motorcycle", "yes");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollMotorcycle.YES, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollMotorcycleNo() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:motorcycle", "no");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollMotorcycle.NO, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollMotorcycleMissing() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollMotorcycle.MISSING, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }
}
