package cytargetlinker.conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapperException;
import org.bridgedb.IDMapperStack;
import org.bridgedb.Xref;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.SpeciesCodes;
import cytargetlinker.conversion.utils.Utils;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesIn;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;

/**
 * Converts miRNA target text files from TarBase to a XGMML or GML network.
 * 
 * The input files for this script can be downloaded here:
 * 
 * 
 * @author Thomas, Pooja, Naime, Tina
 */
public class TarBase {

	private final static Logger log = Logger.getLogger(TarBase.class.getName());
	static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}
	
	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.TarBase 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbSource = BridgeDb mapping file for the source nodes (gene database)
	 * --bridgeDbTarget = BridgeDb mapping file for the target nodes (miRBase database)
	 * --organism "Mouse" or "Human" 
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		TarBase converter = new TarBase();
		converter.startConversion();
	}
	
	private IDMapperStack gdb;
	private Graph graph;
	private Map<String, Integer> index;
	private Map<String, List<String>> edges;
	private List<String> foundConnections;
	private Integer countEdges = 0;
	private boolean mapping = true;

	public TarBase() {
		foundConnections = new ArrayList<String>();
		index = new HashMap<String, Integer>();
		edges = new HashMap<String, List<String>>();
	}

	public void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			Utils.setUpLogger(log, getLogFile(), true);
			if(pargs.isBridgeDbFiles()) gdb = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			if(gdb == null) {
				mapping = false;
				log.info("no identifier mapping");
			}
			if (pargs.getOrganism() != null) {
				log.info("conversion of TarBase file started ...\n");
				ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
					public Graph buildGraph(File in) throws Exception {
						return importMiRecords(in);
					}
				});
				log.info("conversion of TarBase file finalized ...\n");
			} else {
				log.severe("Please specify organism. --organism Human or Mouse.\n");
			}
		} else {
			System.out.println("Please specify input file!");
		}
	}

	public Graph importMiRecords(File input) throws IOException {

		BufferedReader br = new BufferedReader(new FileReader(input));
		String[] header = br.readLine().split("\t");
		index = new HashMap<String, Integer>();
		for (int i = 0; i < header.length; i++) {
			index.put(header[i], i);
		}
		graph = new Graph();
		setNetworkAttributes(input);

		// Load each line into memory
		String line = null;
		List<String[]> rows = new ArrayList<String[]>();
		while ((line = br.readLine()) != null) {
			String[] str = line.split("\t", header.length);
			if (str[index.get("Organism")].equals(pargs.getOrganism())) {
				rows.add(line.split("\t", header.length));
			}
		}

		edges = new HashMap<String, List<String>>();
		for (String[] r : rows) {

			String gene = r[index.get("Ensembl")];
			String mirna = getMiRNA(r[index.get("miRNA")], r[index.get("Organism")]);

			if (gene != null && mirna != null) {
				createGeneNode(gene, r);
				createMiRNANode(mirna, r);
				addEdge(gene, mirna, r);
			}
		}

		log.info(foundConnections.size() + " interactions have been found.\n");

		return graph;
	}
	
	private void createMiRNANode(String id, String [] row) {
		String identifiers = "[" + id;
		String type = "microRNA";
		String mimat = "";
		if(mapping) {
			Xref xrefIn = new Xref(id, DataSource.getBySystemCode("Mb"));
			try {
				Set<Xref> res = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mbm"));
				if(!res.isEmpty()) {
					String str = res.iterator().next().getId();
					if(str.contains("MIMAT")) mimat = str;
					identifiers = identifiers + "," + str;
				}
				Set<Xref> res2 = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mb"));
				if(!res2.isEmpty()) {
					String str = res2.iterator().next().getId();
					if(!str.equals(id)) identifiers = identifiers + "," + str;
				}
			} catch (IDMapperException ex) {
				// could not be mapped
			}
			
		}
		identifiers = identifiers + "]";
		
		Node node = graph.addNode(id);
		node.appendAttribute("identifiers", identifiers);
		node.appendAttribute("label", id);
		node.appendAttribute("name", id);
		node.appendAttribute("miRBaseAccession", mimat);
		node.appendAttribute("biologicalType", type);
		node.appendAttribute("organism", row[index.get("Organism")]);
		node.appendAttribute("Pathology_or_Event", row[index.get("Pathology_or_Event")]);
		node.appendAttribute("miRNA_seq", row[index.get("miRNA_seq")]);
		node.appendAttribute("Seq_location", row[index.get("Seq_location")]);
		node.appendAttribute("Target_seq", row[index.get("Target_seq")]);
	}
	
	private void createGeneNode(String id, String [] row) {
		String geneName = row[index.get("Gene")];
		String geneId = id;
		String type = "gene";
		String organism = row[index.get("Organism")];
		
		if(graph.getNode(geneId) == null) {
			String entrez = "";
			String identifiers = "[" + geneId;
			
			if(mapping) {
				Xref xrefIn = new Xref(geneId, DataSource.getBySystemCode("En"));
				Set<Xref> ens;
				try {
					ens = gdb.mapID(xrefIn, DataSource.getBySystemCode("L"));
					if(!ens.isEmpty()) {
						entrez = ens.iterator().next().getId();
						identifiers = identifiers + "," + entrez;
					}
				} catch (IDMapperException e) {
					// could not be mapped
				}
				
			}
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(geneId);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", geneName);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			node.appendAttribute("ensemblID", id);
			node.appendAttribute("entrezGeneID", entrez);
			node.appendAttribute("HGNC_Symbol", row[index.get("HGNC_Symbol")]);
			node.appendAttribute("Isoform", row[index.get("Isoform")]);
			node.appendAttribute("Chromosomal_location", row[index.get("Chr_loc")]);
			node.appendAttribute("HGNC_ID", row[index.get("HGNC_ID")]);
			node.appendAttribute("SwissProt", row[index.get("SwissProt")]);
		}
	}

	private void setNetworkAttributes(File input) {
		graph.setTitle("TarBase");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "TarBase");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.TYPE.getName(), "miRNA targets");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "Ensembl");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "miRBase");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "mirna");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
	}

	private static String getMiRNA(String id, String species) {
		String code = SpeciesCodes.getCode(species);
		if (code != null) return code + "-" + id;
		return id;
	}

	private void addEdge(String gene, String mirna, String[] r) {

		if (edges.containsKey(gene)) {
			if (!edges.get(gene).contains(mirna)) {
				Edge e = graph.addEdge("" + countEdges, graph.getNode(mirna), graph.getNode(gene));
				setEdgeAttributes(e, r);
				
				edges.get(gene).add(mirna);
				foundConnections.add(gene + "\t" + mirna);
				countEdges++;
			}
		} else {
			Edge e = graph.addEdge("" + countEdges, graph.getNode(mirna), graph.getNode(gene));
			setEdgeAttributes(e, r);

			List<String> list = new ArrayList<String>();
			list.add(mirna);
			foundConnections.add(gene + "\t" + mirna);
			edges.put(gene, list);
			countEdges++;
		}
	}
	
	private void setEdgeAttributes(Edge e, String [] r) {
		e.setAttribute("Data_Type", r[index.get("Data_Type")]);
		e.setAttribute("Support_Type", r[index.get("Support_Type")]);
		e.setAttribute("Organism", r[index.get("Organism")]);
		e.setAttribute("MRE", r[index.get("MRE")]);
		e.setAttribute("S_S_S", r[index.get("S_S_S")]);
		e.setAttribute("I_S", r[index.get("I_S")]);
		e.setAttribute("D_S", r[index.get("D_S")]);
		e.setAttribute("Validation", r[index.get("Validation")]);
		e.setAttribute("Paper", r[index.get("Paper")]);
		e.setAttribute("PMID", r[index.get("PMID")]);
		e.setAttribute("Bibliographic_Notes", r[index.get("Bibliographic_Notes")]);
		e.setAttribute("Cell_Line_Used", r[index.get("Cell_Line_Used")]);
	}
	
	private File getLogFile() {
		if(pargs.isLogFile()) {
			return pargs.getLogFile();
		} else {
			if(pargs.isInput()) {
				return new File(pargs.getInput() + ".log");
			} else {
				return null;
			}
		}
	}
}
