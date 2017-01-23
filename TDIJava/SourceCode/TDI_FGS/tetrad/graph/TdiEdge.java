package tetrad.graph;

///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////


import tetrad.util.TetradSerializable;

import java.io.IOException;
import java.io.ObjectInputStream;

/**
 * Represents an edge node1 *-# node2 where * and # are endpoints of type
 * Endpoint--that is, Endpoint.TAIL, Endpoint.ARROW, or Endpoint.CIRCLE.
 * <p>
 * Note that because speed is of the essence, and Edge cannot be compared
 * to an object of any other type; this will throw an exception.
 *
 * @author Joseph Ramsey
 */
public final class TdiEdge extends Edge {
  
    
    
 //=========================CONSTRUCTORS============================//
   public TdiEdge(Node node1, Node node2, Endpoint endpoint1,
                Endpoint endpoint2) {
       super(node1, node2, endpoint1, endpoint2);
   }
   public TdiEdge(TdiEdge edge) {
        super(edge);
   }
   
  
   
        
}




