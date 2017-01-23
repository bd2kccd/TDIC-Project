/*
 * Copyright (C) 2016 University of Pittsburgh.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package tetrad.io;

import tetrad.data.DataSet;
import java.io.IOException;
import java.util.Set;

/**
 * Interface for reading in data.
 *
 * Feb 29, 2016 1:32:03 PM
 *
 * @author Kevin V. Bui (kvb2@pitt.edu)
 */
public interface DataReader {

    /**
     * Read in dataset.
     *
     * @return
     * @throws IOException whenever unable to read file
     */
    public DataSet readInData() throws IOException;

    /**
     * Read in dataset. Excludes any variables from the given set.
     *
     * @param excludedVariables set of variables to exclude
     * @return
     * @throws IOException whenever unable to read file
     */
    public DataSet readInData(Set<String> excludedVariables) throws IOException;

}
