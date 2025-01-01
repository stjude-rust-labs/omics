use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;

pub struct Corpus {
    hm: Arc<Mutex<HashMap<String, usize>>>,

    lookup: Arc<Vec<String>>,
}

impl Corpus {
    pub fn intern(&self, value: &str) -> usize {
        let mut hm = self.hm.lock().unwrap();

        if let Some(entry) = hm.get(value) {
            return *entry;
        }

        let current = hm.len();
        hm.insert(value.to_owned(), current);
        return current;
    }

    pub fn resolve(&self, id: usize) -> Option<String> {
        self.lookup.get(id).cloned()
    }
}
