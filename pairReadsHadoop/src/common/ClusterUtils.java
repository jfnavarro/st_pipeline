// Copyright (C) 2011-2012 CRS4.
//
// This file is part of Seal.
//
// Seal is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Seal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with Seal.  If not, see <http://www.gnu.org/licenses/>.

package common;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.mapred.ClusterStatus;
import org.apache.hadoop.mapred.JobClient;
import org.apache.hadoop.mapred.JobConf;

import java.io.IOException;

public class ClusterUtils
{
	public static final String NUM_RED_TASKS_PROPERTY = "mapred.reduce.tasks"; // XXX: this changes depending on Hadoop version

	public static int getNumberTaskTrackers(Configuration conf) throws IOException
	{
		/* XXX hack to get the ClusterStatus.  To use JobClient it seems I need to
		 * wrap the Configuration with the deprecated JobConf.
		 * Is there a better way?
		 */
		ClusterStatus status = (new JobClient(new JobConf(conf))).getClusterStatus();
		return status.getTaskTrackers();
	}
}
