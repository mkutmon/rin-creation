package cytargetlinker.conversion.data;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Node;

public class MiRNANode {

	private String id;
	private List<String> names;
	private String label;
	private String biologicalType = "microRNA";
	
	private MiRNANode() {
		names = new ArrayList<String>();
	}
	
	public static MiRNANode createMiRNANode(String name, IDMapper mapper, DataSource in) {
		Xref xrefIn = new Xref(name, in);
		try {
			Set<Xref> result = mapper.mapID(xrefIn, DataSource.getBySystemCode("Mbm"));
			if(result.size() == 1) {
				MiRNANode node = new MiRNANode();
				node.setId(result.iterator().next().getId());
				Set<Xref> names = mapper.mapID(xrefIn, DataSource.getBySystemCode("Mb"));
				for(Xref x : names) {
					if(!node.getNames().contains(x.getId())) node.getNames().add(x.getId());
				}
				node.setLabel(name);
				return node;
			}
		} catch(IDMapperException e) {
			return null;
		}
		return null;
	}
	
	public Node getNode(Graph graph) {
		Node node = graph.addNode(getId());
		node.appendAttribute("identifiers", getIdentifiers());
		node.appendAttribute("label", getLabel());
		node.appendAttribute("name", getLabel());
		node.appendAttribute("biologicalType", getBiologicalType());
		return node;
	}
	
	private String getIdentifiers() {
		String ids = "[";
		ids = ids + getId();
		
		for(String name : getNames()) {
			ids = ids + "," + name;
		}
		
		ids = ids + "]";
		return ids;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public List<String> getNames() {
		return names;
	}

	public void setNames(List<String> names) {
		this.names = names;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public String getBiologicalType() {
		return biologicalType;
	}

	public void setBiologicalType(String biologicalType) {
		this.biologicalType = biologicalType;
	}
}
