use std::collections::VecDeque;
use std::num::NonZeroU32;
use std::rc::Rc;

#[derive(Debug, Clone)]
struct Link<T> {
    data: T,
    prev: Option<Rc<Link<T>>>,
}

impl<T> Link<T> {
    pub fn new(data: T) -> Rc<Self> {
        Self { data, prev: None }.into()
    }

    pub fn link(self: &Rc<Self>, data: T) -> Rc<Self> {
        Rc::from(Link {
            data,
            prev: Some(self.clone()),
        })
    }
}

#[derive(Default, Debug, Clone)]
struct TrieNode {
    a: Option<NonZeroU32>,
    c: Option<NonZeroU32>,
    g: Option<NonZeroU32>,
    t: Option<NonZeroU32>,
}

#[derive(Default, Debug, Clone)]
pub struct KmerTrie {
    storage: Vec<TrieNode>,
    k: usize,
}

impl KmerTrie {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            storage: vec![TrieNode::default()],
        }
    }

    pub fn insert(&mut self, word: &str) {
        assert_eq!(word.len(), self.k);

        let mut current_node = 0;
        for c in word.chars() {
            current_node = {
                let mut node = std::mem::take(&mut self.storage[current_node as usize]);
                let x = match c {
                    'A' => node
                        .a
                        .get_or_insert_with(|| {
                            let n =
                                NonZeroU32::new(self.storage.len().try_into().unwrap()).unwrap();
                            self.storage.push(Default::default());
                            n
                        })
                        .get(),
                    'C' => node
                        .c
                        .get_or_insert_with(|| {
                            let n =
                                NonZeroU32::new(self.storage.len().try_into().unwrap()).unwrap();
                            self.storage.push(Default::default());
                            n
                        })
                        .get(),
                    'G' => node
                        .g
                        .get_or_insert_with(|| {
                            let n =
                                NonZeroU32::new(self.storage.len().try_into().unwrap()).unwrap();
                            self.storage.push(Default::default());
                            n
                        })
                        .get(),
                    'T' => node
                        .t
                        .get_or_insert_with(|| {
                            let n =
                                NonZeroU32::new(self.storage.len().try_into().unwrap()).unwrap();
                            self.storage.push(Default::default());
                            n
                        })
                        .get(),
                    _ => unreachable!(),
                };
                self.storage[current_node as usize] = node;
                x
            };
        }
    }

    pub fn contains(&self, word: &str) -> bool {
        assert_eq!(word.len(), self.k);

        let mut current_node = 0;
        for c in word.chars() {
            let option = match c {
                'A' => &self.storage[current_node as usize].a,
                'C' => &self.storage[current_node as usize].c,
                'G' => &self.storage[current_node as usize].g,
                'T' => &self.storage[current_node as usize].t,
                _ => unreachable!(),
            };
            match option {
                Some(node) => current_node = node.get(),
                None => return false,
            }
        }

        true
    }

    pub fn fuzzy_search(&self, word: &str, max_mismatches: usize) -> Vec<(String, usize)> {
        if word.len() < self.k {
            panic!("{} < {}", word.len(), self.k);
        }

        struct State {
            node: u32,
            mismatches: usize,
            position: usize,
            prefix: Rc<Link<u8>>,
        }

        let mut results = Vec::new();
        let mut bytes = word.as_bytes();
        let mut stack = VecDeque::new();
        stack.push_back(State {
            node: 0,
            mismatches: 0,
            position: 0,
            prefix: Link::new(0),
        });

        while let Some(state) = stack.pop_back() {
            if state.position == self.k {
                let mut s = Vec::new();
                let mut n = &state.prefix;

                while let Some(x) = &n.prev {
                    s.push(n.data);
                    n = x;
                }

                s.reverse();
                results.push((String::from_utf8(s).unwrap(), state.mismatches));
                // unimplemented!()
            } else {
                let c = bytes[state.position];
                if let Some(node) = &self.storage[state.node as usize].a {
                    if c == b'A' || state.mismatches < max_mismatches {
                        stack.push_back(State {
                            node: node.get(),
                            mismatches: state.mismatches + (c != b'A' && c != b'N') as usize,
                            position: state.position + 1,
                            prefix: state.prefix.link(b'A'),
                        })
                    }
                }
                if let Some(node) = &self.storage[state.node as usize].c {
                    if c == b'C' || state.mismatches < max_mismatches {
                        stack.push_back(State {
                            node: node.get(),
                            mismatches: state.mismatches + (c != b'C' && c != b'N') as usize,
                            position: state.position + 1,
                            prefix: state.prefix.link(b'C'),
                        })
                    }
                }
                if let Some(node) = &self.storage[state.node as usize].g {
                    if c == b'G' || state.mismatches < max_mismatches {
                        stack.push_back(State {
                            node: node.get(),
                            mismatches: state.mismatches + (c != b'G' && c != b'N') as usize,
                            position: state.position + 1,
                            prefix: state.prefix.link(b'G'),
                        })
                    }
                }
                if let Some(node) = &self.storage[state.node as usize].t {
                    if c == b'T' || state.mismatches < max_mismatches {
                        stack.push_back(State {
                            node: node.get(),
                            mismatches: state.mismatches + (c != b'T' && c != b'N') as usize,
                            position: state.position + 1,
                            prefix: state.prefix.link(b'T'),
                        })
                    }
                }
            }
        }

        results
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_contains() {
        let mut trie = KmerTrie::new(4);
        assert!(!trie.contains("ATTA"));
        assert!(!trie.contains("ATTT"));

        trie.insert("ATTA");
        trie.insert("ATGC");
        assert!(!trie.contains("ATTC"));
        assert!(!trie.contains("CGTC"));
        assert!(trie.contains("ATTA"));
        assert!(trie.contains("ATGC"));
    }

    #[test]
    fn test_fuzzy_search() {
        let mut trie = KmerTrie::new(4);
        trie.insert("ATTA");
        trie.insert("ATTT");
        trie.insert("ATGC");

        let mut results = trie.fuzzy_search("ATTA", 1);
        println!("{:?}", results);
        results.sort();

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].0, "ATTA");
        assert_eq!(results[1].0, "ATTT");
        assert_eq!(results[0].1, 0);
        assert_eq!(results[1].1, 1);
    }
}
