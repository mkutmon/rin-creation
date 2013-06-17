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
 * Converts miRNA target text files from MiRecords to a XGMML or GML network.
 * 
 * The input files for this script can be downloaded here:
 * http://mirecords.biolead.org/download_data.php?v=2
 * 
 * @author Thomas, Pooja, Naime, Tina
 */
public class MiRecords {

	private final static Logger log = Logger.getLogger(MiRecords.class.getName());
	static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {}

	private String dbName = "miRecords";
	private String dbVersion;
	
	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.MiRecords 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbFile = BridgeDb mapping files for source and target nodes
	 * --organism "Mus musculus" or "Homo sapiens" 
	 * -v = database version
	 */
	public static void main(String argv[]) throws Exception {

		pargs = ArgsParser.parse(argv, Args.class);

		MiRecords converter = new MiRecords();
		converter.startConversion(pargs);
	}

	private IDMapperStack gdb;
	private Graph graph;
	private Map<String, Integer> index;
	private Map<String, List<String>> edges;

	private List<String> foundConnections;
	private Integer countEdges = 0;
	private boolean mapping = true;
	private Integer countGenes = 0;
	private Integer countMiRNAs = 0;

	public MiRecords() {
		index = new HashMap<String, Integer>();
		edges = new HashMap<String, List<String>>();
		foundConnections = new ArrayList<String>();
	}

	public void startConversion(Args pargs) throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			if(pargs.isDatabaseVersion()) {
				dbVersion = pargs.getDatabaseVersion();
			}
			Utils.setUpLogger(log, getLogFile(), true);
			if(pargs.isBridgeDbFiles()) gdb = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			if(gdb == null) {
				mapping = false;
				log.info("no identifier mapping");
			}
			if (pargs.getOrganism() != null) {
				log.info("conversion of MiRecords file started ...\n");
				ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
					public Graph buildGraph(File in) throws Exception {
						return importMiRecords(in);
					}
				});
	
				log.info("miRecords network file created.\n\n");
			} else {
				log.severe("Please specify organism. --organism Homo sapiens or Mus musculus.\n");
			}
		} else {
			System.out.println("Please specify input file!");
		}
	}

	public Graph importMiRecords(File input) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(input));
		String[] header = br.readLine().split("\t");

		for (int i = 0; i < header.length; i++)
			index.put(header[i], i);

		graph = new Graph();
		setNetworkAttributes(input);

		// Load each line into memory
		String line = null;
		List<String[]> rows = new ArrayList<String[]>();
		while ((line = br.readLine()) != null) {
			line = removeInvalidXMLCharacters(line);
			String[] str = line.split("\t", header.length);
			if (str[index.get("Target gene_species_scientific")].equals(pargs.getOrganism())) {
				rows.add(line.split("\t", header.length));
			}
		}

		for (String[] row : rows) {
			String geneNode = createGeneNode(row);
			String miRNANode = createMiRNANode(row);
			addEdge(geneNode, miRNANode, row);
		}

		log.info(foundConnections.size() + " interactions have been found.\n" + countGenes + " gene nodes.\n" + countMiRNAs + " miRNA nodes.\n");
		cleanUp();
		return graph;
	}
	
	private String createGeneNode(String[] r) {
		String geneName = r[index.get("Target gene_name")];
		String geneId = r[index.get("Target gene_Refseq_acc")];
		String type = "gene";
		String organism = r[index.get("Target gene_species_scientific")];
		
		if(geneId.contains(".")) {
			int pos = geneId.indexOf(".");
			geneId = geneId.substring(0, pos);
		}
		
		if(graph.getNode(geneId) == null) {
			String ensembl = "";
			String entrez = "";
			String identifiers = "[" + geneId;
			
			if(mapping) {
				Xref xrefIn = new Xref(geneId, DataSource.getBySystemCode("Q"));
				try {
					Set<Xref>ens = gdb.mapID(xrefIn, DataSource.getBySystemCode("En"));
					if(!ens.isEmpty()) {
						ensembl = ens.iterator().next().getId();
						identifiers = identifiers + "," + ensembl;
					}
					Set<Xref>entrezRes = gdb.mapID(xrefIn, DataSource.getBySystemCode("L"));
					if(!entrezRes.isEmpty()) {
						entrez = entrezRes.iterator().next().getId();
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
			node.appendAttribute("name", geneName);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			node.appendAttribute("ensemblID", ensembl);
			node.appendAttribute("entrezGeneID", entrez);
			node.appendAttribute("refseq", geneId);
			countGenes++;
		}
		
		return geneId;
	}
	
	private String createMiRNANode(String[] r) {
		String miRNA = processMirna(r[index.get("miRNA_mature_ID")], r[index.get("miRNA_species")]);
		
		String type = "microRNA";
		String organism = r[index.get("miRNA_species")];
		
		if(graph.getNode(miRNA) == null) {
			String mimat = "";
			String identifiers = "[" + miRNA;
			Xref xrefIn = new Xref(miRNA, DataSource.getBySystemCode("Mb"));
			
			if(mapping) {
				try {
					Set<Xref> result = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mb"));
					Set<Xref> result2 = gdb.mapID(xrefIn, DataSource.getBySystemCode("Mbm"));
					
					List<String> list = new ArrayList<String>();
					list.add(miRNA);
					for(Xref x : result) {
						if(!list.contains(x.getId())) {
							list.add(x.getId());
							identifiers = identifiers + "," + x.getId();
						}
					}
					for(Xref x : result2) {
						if(!list.contains(x.getId())) {
							list.add(x.getId());
							if(x.getId().startsWith("MIMAT")) {
								mimat = x.getId();
							}
							identifiers = identifiers + "," + x.getId();
						}
					}
				} catch (IDMapperException e) {
					// could not be mapped
				}
			}
			identifiers = identifiers + "]";
			
			Node node = graph.addNode(miRNA);
			node.appendAttribute("identifiers", identifiers);
			node.appendAttribute("label", miRNA);
			node.appendAttribute("name", miRNA);
			node.appendAttribute("miRBaseAccession", mimat);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			countMiRNAs++;
		}
		
		return miRNA;
	}
	

	private void cleanUp() {
		edges.clear();
		index.clear();
		foundConnections.clear();
	}

	private void setNetworkAttributes(File input) {
		graph.setTitle("MiRecords");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), dbName + " " + dbVersion + " (" + pargs.getOrganism() + ")");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "MTI");
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), input.getName());
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), "RefSeq");
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "miRBase");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "mirna");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
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

	private void setEdgeAttributes(Edge e, String[] r) {
		e.setAttribute("miRNA_regulation", r[index.get("miRNA_regulation")]);
		e.setAttribute("Reporter_target_gene_region", r[index.get("Reporter_target gene/region")]);
		e.setAttribute("Reporter_link_element", r[index.get("Reporter link element")]);
		e.setAttribute("Test_method_inter", r[index.get("Test_method_inter")]);
		e.setAttribute("Target_gene_mRNA_level", r[index.get("Target gene mRNA_level")]);
		e.setAttribute("Description", r[index.get("Original description")]);
		e.setAttribute("Mutation_target_region", r[index.get("Mutation_target region")]);
		e.setAttribute("Post_mutation_method", r[index.get("Post mutation_method")]);
		e.setAttribute("Original_description_mutation_region", r[index.get("Original description_mutation_region")]);
		e.setAttribute("Target_site_position", r[index.get("Target site_position")]);
		e.setAttribute("miRNA_regulation_site", r[index.get("miRNA_regulation_site")]);
		e.setAttribute("Reporter_target_site", r[index.get("Reporter_target site")]);
		e.setAttribute("Reporter_link_element", r[index.get("Reporter link element")]);
		e.setAttribute("Test_method_inter_site", r[index.get("Test_method_inter_site")]);
		e.setAttribute("Original_description_inter_site", r[index.get("Original description_inter_site")]);
		e.setAttribute("Mutation_target_site", r[index.get("Mutation_target site")]);
		e.setAttribute("Post_mutation_method_site", r[index.get("Post mutation_method_site")]);
		e.setAttribute("Original_description_mutation_site", r[index.get("Original description_mutation_site")]);
		e.setAttribute("Additional_note", r[index.get("Additional note")]);
		e.setAttribute("Pubmed_id", r[index.get("Pubmed_id")]);
		e.setAttribute("interactionType", "MTI");
		e.setAttribute("datasource", dbName + " " + dbVersion);
	}

	private String processMirna(String mirna, String species) {
		mirna = mirna.replace("[", "");
		mirna = mirna.replace("]", "");

		if (mirna.startsWith("miR")) {
			// Add species code
			String code = SpeciesCodes.getCode(species);
			if (code != null)
				mirna = code + "-" + mirna;
		}
		return mirna;
	}

	/**
	 * Removes all invalid Unicode characters that are not suitable to be used
	 * either in markup or text inside XML Documents.
	 * 
	 * Based on these recommendations
	 * http://www.w3.org/TR/2000/REC-xml-20001006#NT-Char
	 * http://cse-mjmcl.cse.bris.ac.uk/blog/2007/02/14/1171465494443.html
	 * 
	 * @param s: The resultant String stripped of the offending characters!
	 * @return
	 */
	public String removeInvalidXMLCharacters(String s) {
		StringBuilder out = new StringBuilder();

		int codePoint;
		int i = 0;

		while (i < s.length()) {
			// This is the unicode code of the character.
			codePoint = s.codePointAt(i);
			if ((codePoint == 0x9) || (codePoint == 0xA) || (codePoint == 0xD)
					|| ((codePoint >= 0x20) && (codePoint <= 0xD7FF))
					|| ((codePoint >= 0xE000) && (codePoint <= 0xFFFD))
					|| ((codePoint >= 0x10000) && (codePoint <= 0x10FFFF))) {
				out.append(Character.toChars(codePoint));
			}
			i += Character.charCount(codePoint);
		}
		return out.toString();
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
