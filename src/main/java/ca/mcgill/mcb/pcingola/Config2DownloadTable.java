package ca.mcgill.mcb.pcingola;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;

/**
 * Create an HTML 'download' table based on the config file
 * Also creates a list of genome for Galaxy menu
 * 
 * @author pablocingolani
 */
public class Config2DownloadTable extends SnpEff {

	public static final String DARK_ROW = "bgcolor=#CCCCCC";
	public static final String LIGHT_ROW = "bgcolor=#EEEEEE";

	public static final String HTTP_PROTOCOL = "http://";
	public static final String FTP_PROTOCOL = "ftp://";

	boolean galaxy = false;
	Config config;
	HashSet<String> names;
	ArrayList<String> namesSorted;
	ArrayList<String> genVerSorted;

	public static void main(String[] args) {
		Config2DownloadTable conf2down = new Config2DownloadTable();
		conf2down.parseArgs(args);
		conf2down.run();
	}

	public Config2DownloadTable() {
		// Read config (it doesn't matter which genome)
		config = new Config("hg19");

		// Get all genome names and sort them
		names = new HashSet<String>();
		for (String genVer : config)
			names.add(config.getName(genVer));

		namesSorted = new ArrayList<String>();
		namesSorted.addAll(names);
		Collections.sort(namesSorted);

		// Sort genome versions
		genVerSorted = new ArrayList<String>();
		for (String genVer : config)
			genVerSorted.add(genVer);
		Collections.sort(genVerSorted, Collections.reverseOrder());
	}

	/**
	 * Create download table
	 */
	void downloadTable() {
		// Create an HTML table
		boolean dark = false;
		String bg = "";

		System.out.println("\t<table> <tr " + DARK_ROW + "> <td> <b> Genome </b> </td>  <td> <b> Version </b> </td>  <td> <b> Reference </b> </td> </tr>");
		for (String name : namesSorted) {

			// Color
			if (dark) bg = DARK_ROW;
			else bg = LIGHT_ROW;
			dark = !dark;

			boolean showName = true;
			for (String genVer : genVerSorted) {
				String n = config.getName(genVer);
				// In this group?
				if (name.equals(n)) {
					System.out.println("\t\t<tr " + bg + ">");

					// Show name
					String name2show = showName ? name.replace('_', ' ') : "&nbsp;";
					System.out.println("\t\t\t<td> " + name2show + " </td>");
					showName = false;

					// Download link
					String url = "http://sourceforge.net/projects/snpeff/files/databases/v" + SnpEff.VERSION_MAJOR + "/snpEff_v" + SnpEff.VERSION_MAJOR + "_" + genVer + ".zip";
					System.out.println("\t\t\t<td> <a class=\"body\" href=\"" + url + "\"> " + genVer + " </a> </td>");

					// Reference
					String ref = config.getReference(genVer);
					String link = "";
					if (ref != null) {
						if (ref.indexOf(',') > 0) ref = ref.substring(0, ref.indexOf(',')); // Many references? Use the first one
						link = ref;

						// Remove everything after slash
						int idx = ref.indexOf('/', HTTP_PROTOCOL.length());
						if (idx > 0) ref = ref.substring(0, idx);

						// Remove protocol
						if (ref.startsWith(HTTP_PROTOCOL)) ref = ref.substring(HTTP_PROTOCOL.length());
						if (ref.startsWith(FTP_PROTOCOL)) ref = ref.substring(FTP_PROTOCOL.length());

					} else ref = "";
					System.out.println("\t\t\t<td> <a class=\"body\" href=\"" + link + "\">" + ref + "</a> </td>");

					System.out.println("\t\t</tr>");
				}
			}
		}
		System.out.println("\t</table>");
	}

	/**
	 * Galaxy config genome list
	 */
	void galaxyConfig() {
		System.out.println("\t<param name=\"genomeVersion\" type=\"select\" label=\"Genome\">");

		for (String name : namesSorted) {
			for (String genVer : genVerSorted) {
				String n = config.getName(genVer);

				// In this group?
				if (name.equals(n)) {
					System.out.println("\t\t<option value=\"" + genVer + "\">" + name.replace('_', ' ') + " : " + genVer + "</option>");
				}
			}
		}

		System.out.println("\t</param>");
	}

	@Override
	public void parseArgs(String[] args) {
		if (args.length != 1) usage(null);
		galaxy = args[0].equals("galaxy");
	}

	@Override
	public boolean run() {
		if (galaxy) galaxyConfig();
		else downloadTable();

		return true;
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("Usage: snpEff cfg2table [download | galaxy]");
		System.exit(-1);
	}
}
