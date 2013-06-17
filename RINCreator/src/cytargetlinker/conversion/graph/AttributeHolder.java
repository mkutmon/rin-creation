package cytargetlinker.conversion.graph;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * 
 * @author Thomas
 *
 */
public class AttributeHolder {
	Map<String, Object> attributes = new HashMap<String, Object>();

	public void setAttribute(String name, String value) {
		attributes.put(name, value);
	}

	public Object getAttribute(String name) {
		return attributes.get(name); 
	}

	public Set<String> getAttributeNames() {
		return attributes.keySet();
	}

	public void appendAttribute(String name, String value) {
		appendAttribute(name, value, "; ");
	}
	
	public void appendAttribute(String name, String value, String sep) {
		String curr = attributes.get(name) == null ? "" : attributes.get(name).toString();
		if("".equals(curr)) curr = value;
		else if(!curr.startsWith(value) && !curr.contains(sep + value)) {
			curr += sep + value;
		}
		attributes.put(name, curr);
	}
}
