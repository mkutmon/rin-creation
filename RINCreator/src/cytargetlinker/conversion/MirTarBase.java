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

import cytargetlinker.conversion.data.MiRNANode;
import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesIn;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.Utils;

/**
 * Converts miRNA target text files from MiRTarBase to a XGMML or GML network.
 * 
 * The input files for this script can be downloaded here:
 * http://mirtarbase.mbc.nctu.edu.tw/php/download.php (under
 * "Catatolog by Species").
 * 
 * @author Thomas, Pooja, Naime, Tina
 */
public class MirTarBase {

	private String dbName = "miRTarBase";
	private String dbVersion;
			
	private final static Logger log = Logger.getLogger(MirTarBase.class.getName());
	static Args pargs;

	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}

	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.MirTarBase 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbFile = BridgeDb mapping files for source and target nodes
	 * --organism like listed in the interaction file "Mus musculus" or "Homo sapiens" 
	 * -v database version (2.5 / 3.5)
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		MirTarBase converter = new MirTarBase();
		converter.startConversion();
	}

	private Graph graph;
	private Map<String, Integer> index;
	private Map<String, List<String>> edges;
	private IDMapperStack gdb;

	private List<String> foundConnections;
	private Integer countEdges = 0;
	private boolean mapping = true;
	private Integer countGenes = 0;
	private Integer countMiRNAs = 0;

	public MirTarBase() {
		foundConnections = new ArrayList<String>();
		index = new HashMap<String, Integer>();
		edges = new HashMap<String, List<String>>();
	}

	public void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			if(pargs.isDatabaseVersion()) {
				dbVersion = pargs.getDatabaseVersion();
			}
			Utils.setUpLogger(log, getLogFile(), false);
			log.info("conversion of miRTarBase file started ...\n");
			if(pargs.isBridgeDbFiles()) gdb = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			
			if (gdb == null) {
				mapping = false;
				log.info("no identifier mapping");
			}
			if (pargs.getOrganism() != null) {
				ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
					public Graph buildGraph(File in) throws Exception {
						return createNetwork(in);
					}
				});
			} else {
				log.severe("Please specify organism. --organism Homo sapiens or Mus musculus.\n");
			}
			log.info("miRTarBase network file created.\n\n");
		} else {
			System.out.println("Please specify input file!");
		}
	}

	public Graph createNetwork(File input) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(input));

		String[] header = br.readLine().split("\t");
		index = new HashMap<String, Integer>();
		for (int i = 0; i < header.length; i++) {
			index.put(header[i], i);
		}

		graph = new Graph();
		setNetworkAttributes(input);

		String line = null;
		List<String[]> rows = new ArrayList<String[]>();
		while ((line = br.readLine()) != null) {
			String[] str = line.split("\t");
			if (str[index.get("Species (Target Gene)")].equals(pargs .getOrganism())) {
				rows.add(line.split("\t", header.length));
			}
		}
		
		edges = new HashMap<String, List<String>>();
		for (String[] r : rows) {
			String geneNode = createGeneNode(r);
			String miRNANode = createMiRNANode(r);
			if(geneNode != null && miRNANode != null) addEdge(geneNode, miRNANode, r);	
		}

		log.info(foundConnections.size() + " interactions have been found.\n" + countGenes + " gene nodes.\n" + countMiRNAs + " miRNA nodes.\n");
		
		return graph;
	}

	private String createGeneNode(String[] r) {
		String geneName = r[index.get("Target Gene")];
		String geneId = r[index.get("Target Gene (Entrez ID)")];
		String type = "gene";
		String organism = r[index.get("Species (Target Gene)")];
		
		if(graph.getNode(geneId) == null) {
			String identifiers = "[" + geneId;
			String ensembl = "";
			if(mapping) {
				Xref xrefIn = new Xref(geneId, DataSource.getBySystemCode("L"));
				Set<Xref> ens;
				try {
					ens = gdb.mapID(xrefIn, DataSource.getBySystemCode("En"));
					ens.addAll(gdb.mapID(xrefIn, DataSource.getBySystemCode("S")));
					for(Xref x : ens) {
						if(ensembl.equals("") && x.getDataSource().getSystemCode().equals("En")) ensembl = x.getId();
						String id = x.getId();
						identifiers = identifiers + "," + id;
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
			node.appendAttribute("ensemblID", ensembl);
			node.appendAttribute("entrezGeneID", geneId);
			countGenes++;
		}
		
		return geneId;
	}
	
	private String createMiRNANode(String[] r) {
		String miRNA = r[index.get("miRNA")];
		String organism = r[index.get("Species (miRNA)")];

		MiRNANode mirnaNode = MiRNANode.createMiRNANode(miRNA, gdb, DataSource.getBySystemCode("Mb"));
		if(mirnaNode != null) {
			Node node = mirnaNode.getNode(graph);
			node.appendAttribute("miRBaseAccession", mirnaNode.getId());
			node.appendAttribute("organism", organism);
			countMiRNAs++;			
			return mirnaNode.getId();
		} else {
			log.warning("miRNA node for " + miRNA + " not found! Skipping this miRNA.");
			return null;
		}
		
	}

	private void addEdge(String gene, String mirna, String[] r) {
		if (edges.containsKey(gene)) {
			if (!edges.get(gene).contains(mirna)) {
				Edge e = graph.addEdge("" + countEdges, graph.getNode(mirna), graph.getNode(gene));
				e.setAttribute("experiments", r[index.get("Experiments")]);
				e.setAttribute("supportType", r[index.get("Support Type")]);
				e.setAttribute("referenceID", r[index.get("References (PMID)")]);
				e.setAttribute("interactionType", "MTI");
				e.setAttribute("datasource", dbName + " " + dbVersion);
				e.setAttribute("miRTarBaseID", r[index.get("miRTarBase ID")]);
				edges.get(gene).add(mirna);
				foundConnections.add(gene + "\t" + mirna);
				countEdges++;
			}
		} else {
			Edge e = graph.addEdge("" + countEdges, graph.getNode(mirna), graph.getNode(gene));
			e.setAttribute("experiments", r[index.get("Experiments")]);
			e.setAttribute("supportType", r[index.get("Support Type")]);
			e.setAttribute("referenceID", r[index.get("References (PMID)")]);
			e.setAttribute("interactionType", "MTI");
			e.setAttribute("datasource", dbName + " " + dbVersion);
			e.setAttribute("miRTarBaseID", r[index.get("miRTarBase ID")]);
			
			List<String> list = new ArrayList<String>();
			list.add(mirna);
			foundConnections.add(gene + "\t" + mirna);
			edges.put(gene, list);
			countEdges++;
		}
	}
	
	private void setNetworkAttributes(File input) {
		graph.setTitle("MirTarBase");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), dbName + " " + dbVersion + " (" + pargs.getOrganism() + ")");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "MTIs");
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "Entrez Gene");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "miRBase");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "mirna");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
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
