package cytargetlinker.conversion;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.IDMapperStack;
import org.bridgedb.Xref;
import org.pathvisio.core.model.Pathway;

import uk.co.flamingpenguin.jewel.cli.Option;
import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.ArgsParser;
import cytargetlinker.conversion.utils.CommonAttributes;
import cytargetlinker.conversion.utils.Utils;
import cytargetlinker.conversion.utils.ArgsParser.AFilesAttributes;
import cytargetlinker.conversion.utils.ArgsParser.AFilesIn;
import cytargetlinker.conversion.utils.ArgsParser.AFilesOut;
import cytargetlinker.conversion.utils.ArgsParser.AHelp;
import cytargetlinker.conversion.utils.ArgsParser.GraphBuilder;

public class WP2XgmmlConverter {

	private final static Logger log = Logger.getLogger(WP2XgmmlConverter.class.getName());
	private static Args pargs;
	private interface Args extends AHelp, AFilesIn, AFilesOut, AFilesAttributes {
		@Option(longName = "source", description = "Source to use in the network metadata")
		public String getSource();
		public boolean isTitle();
		@Option(longName = "date", description = "Date to use in the network metadata")
		public String getDate();
		public boolean isDate();
		@Option(longName = "mapToCode", description = "Map all genes to specified system")
		public String getOutCode();
		public boolean isOutCode();
	}
	
	/**
	 * MAIN METHOD
	 * USAGE: java -cp conversion.jar cytargetlinker.conversion.WP2XgmmlConverter 
	 * ARGUMENTS: 
	 * -i = input file
	 * -o = output file
	 * --bridgeDbFile = BridgeDb mapping files for source and target nodes
	 * --organism "Homo sapiens" or "Mus musculus" 
	 * @throws Exception 
	 */
	public static void main(String argv[]) throws Exception {
		pargs = ArgsParser.parse(argv, Args.class);
		if(pargs.getOrganism() != null) {
			WP2XgmmlConverter converter = new WP2XgmmlConverter();
			converter.startConversion();
		} else {
			log.severe("Please specify organism (--organism Homo sapiens or Mus musculus).");
		}
	}
	
	public WP2XgmmlConverter() {
		edges = new HashMap<String, List<String>>();
		foundConnections = new ArrayList<String>();
	}
	
	private IDMapperStack gdb;
	private Graph graph;
	private Map<String, List<String>> edges;
	
	private List<String> foundConnections;
	
	private int countEdge = 0;
	private boolean mapping = true;
	
	private void startConversion() throws Exception {
		File in = pargs.getInput();
		if(in != null) {
			Utils.setUpLogger(log, getLogFile(), false);
			if(pargs.isBridgeDbFiles()) gdb = Utils.initIDMapper(pargs.getBridgeDbFiles(), false);
			if(gdb == null) {
				mapping = false;
				log.info("no identifier mapping");
			}
			log.info("Conversion of WikiPathways collection started\n");

			ArgsParser.convertAndWrite(pargs, pargs, new GraphBuilder() {
				public Graph buildGraph(File in) throws Exception {
					return parseWikiPathways(in);
				}
			});
				
			log.info("WikiPathways association network file created.\n\n");
		} else {
			System.out.println("Please specify input file!");
		}
	}
	
	private String getSysCodeIn() {
		if(pargs.isOutCode()) {
			return pargs.getOutCode();
		} else {
			return "En";
		}
	}
	
	private Graph parseWikiPathways(File in) {
		graph = new Graph();
		
		setNetworkAttributes(in);
		parsePathwayCollection(in);
		
		log.info(foundConnections.size() + " interactions have been found.\n" + countGene + " gene and " + countPathway + " pathway nodes.\n");
		return graph;
	}

	private void parsePathwayCollection(File in) {
		if(in.isDirectory()) {
			for (File file : in.listFiles()) {
				try {
					List<String> list = new ArrayList<String>();
					Pathway pathway = new Pathway();
					pathway.setSourceFile(file);
					pathway.readFromXml(file, false);

					String title = pathway.getMappInfo().getMapInfoName();
					
					log.info("Start conversion of pathway " + title + "\n");
					createPathwayNode(pathway.getSourceFile().getName(), title);
					for(Xref xref : pathway.getDataNodeXrefs()) {
						if(xref != null &&  xref.getId() != null && xref.getDataSource() != null && !xref.getId().equals("") && !xref.getDataSource().equals("")) {
							Set<Xref> result = gdb.mapID(xref, DataSource.getBySystemCode(getSysCodeIn()));
							if(result.size() > 0) {
								for(Xref x : result) {
									if(!list.contains(x.getId())) {
										createGeneNode(x, gdb);
										addEdge(x.getId(), pathway.getSourceFile().getName());
										list.add(x.getId());
									}
								}
							}
						}
					}
				} catch (Exception e) {
					log.info("Couldn't convert " + file.getAbsolutePath() + "\n");
				}
			} 
		}
	}
	
	private void addEdge(String gene, String pathway) {
		if(edges.containsKey(gene)) {
			if(!edges.get(gene).contains(pathway)) {
				graph.addEdge("" + countEdge, graph.getNode(pathway),graph.getNode(gene));
				edges.get(gene).add(pathway);
				foundConnections.add(gene + "\t" + pathway);
				countEdge++;
			}
		} else {
			graph.addEdge("" + countEdge, graph.getNode(pathway),graph.getNode(gene));
			List<String> list = new ArrayList<String>();
			list.add(pathway);
			foundConnections.add(gene + "\t" + pathway);
			edges.put(gene, list);
			countEdge++;
		}
	}
	
	private int countGene = 0;
	
	private void createGeneNode(Xref xref, IDMapper mapper) {
		String geneName = xref.getId();
		String geneId = xref.getId();
		String type = "gene";
		String organism = pargs.getOrganism();
		String ensembl = xref.getId();
		
		Map<String, String> attr = new HashMap<String, String>();
		attr.put("label", geneId);
		attr.put("name", geneId);
		attr.put("organism", pargs.getOrganism());
		attr.put("geneid", geneId);
		if(graph.getNode(geneId) == null) {
			String entrez = "";
			String identifiers = "[" + ensembl;
			
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
			node.appendAttribute("name", geneName);
			node.appendAttribute("biologicalType", type);
			node.appendAttribute("organism", organism);
			node.appendAttribute("ensemblID", ensembl);
			node.appendAttribute("entrezGeneID", entrez);
			node.appendAttribute("geneid", geneId);
			countGene++;
		}
	}
	
	private int countPathway = 0;
	
	private Node createPathwayNode(String id, String name) {
		if(graph.getNode(id) == null) {
			String type = "pathway";
			
			Node node = graph.addNode(id);
			node.appendAttribute("identifiers", "[" + id + "]");
			node.appendAttribute("label", name);
			node.appendAttribute("name", name);
			
			node.appendAttribute("biologicalType", type);
			countPathway++;
			return node;
		} else {
			return graph.getNode(id);
		}
	}

	private void setNetworkAttributes(File in) {
		String src = pargs.isSource() ? pargs.getSource() : "WikiPathways";
		String date = pargs.isDate() ? pargs.getDate() : "ND";
		
		graph.setTitle(src + " " + pargs.getOrganism());
		graph.setAttribute(CommonAttributes.DATABASE.getName(), src + " " + pargs.getOrganism());
		graph.setAttribute(CommonAttributes.TYPE.getName(), "pathway annotations");
		graph.setAttribute(CommonAttributes.SOURCE_FILE.getName(), src + " " + pargs.getOrganism() + " " + date);
		graph.setAttribute(CommonAttributes.IDENTIFIERS.getName(), "identifiers");
		graph.setAttribute(CommonAttributes.TARGET_DATASOURCE.getName(), src);
		graph.setAttribute(CommonAttributes.SOURCE_DATASOURCE.getName(), "geneid");
		graph.setAttribute(CommonAttributes.SOURCE_TYPE.getName(), "pathway");
		graph.setAttribute(CommonAttributes.TARGET_TYPE.getName(), "gene");
		graph.setAttribute("Download-Date", date);
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
