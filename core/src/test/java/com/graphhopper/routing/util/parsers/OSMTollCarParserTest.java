package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.*;
import com.graphhopper.storage.IntsRef;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class OSMTollCarParserTest {
    private EnumEncodedValue<TollCar> tollEnc;
    private OSMTollCarParser parser;

    @BeforeEach
    public void setUp() {
        tollEnc = new EnumEncodedValue<>(TollCar.KEY, TollCar.class);
        tollEnc.init(new EncodedValue.InitializerConfig());
        parser = new OSMTollCarParser(tollEnc);
    }

    @Test
    public void testTollCarYes() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:car", "yes");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollCar.YES, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollCarNo() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:car", "no");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollCar.NO, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollCarMissing() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollCar.MISSING, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }
}
