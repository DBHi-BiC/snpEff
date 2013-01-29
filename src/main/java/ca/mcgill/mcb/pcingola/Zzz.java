package ca.mcgill.mcb.pcingola;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.logStatsServer.LogStats;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEffCmdDownload;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * Simple test program
 * @author pcingola
 */
public class Zzz {

	public static void main(String[] args) {
		HashMap<String, String> reportValues = new HashMap<String, String>();
		LogStats logStats = LogStats.report(SnpEff.SOFTWARE_NAME, SnpEff.VERSION_SHORT, SnpEff.VERSION, true, true, args, "", reportValues);

		if (logStats.isNewVersion()) {
			Gpr.debug("New version found: " //
					+ "\n\tNew version  : " + logStats.getLatestVersion() // 
					+ "\n\tRelease date : " + logStats.getLatestReleaseDate() //
					+ "\n\tDownload URL : " + logStats.getLatestUrl() //
			);

			// Invoke download command
			SnpEffCmdDownload download = new SnpEffCmdDownload();
			String argsDownload[] = { "-v", "snpeff" };
			download.parseArgs(argsDownload);
			download.run();
		}

		Gpr.debug("DONE!");
	}
}
