const API = "/api";

function showNotification(message, type = "success") {
  const container = document.getElementById("notifications");
  const notification = document.createElement("div");
  notification.className = `notification ${type}`;
  notification.textContent = message;
  container.appendChild(notification);
  setTimeout(() => notification.remove(), 4000);
}

async function addMolecule() {
  const body = {
    smiles: document.getElementById("smiles").value.trim(),
    name: document.getElementById("name").value.trim() || null,
    id: document.getElementById("id").value.trim() || null,
  };

  try {
    const response = await fetch(`${API}/molecules`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });

    if (response.ok) {
      const data = await response.json();
      showNotification(`✅ Молекула "${data.name || data.smiles}" добавлена`, "success");
      document.getElementById("addForm").reset();
      loadMolecules();
    } else {
      const error = await response.json();
      showNotification(`❌ Ошибка: ${error.detail || "Неизвестная ошибка"}`, "error");
    }
  } catch (err) {
    showNotification(`❌ Ошибка сети: ${err.message}`, "error");
  }
}

async function loadMolecules() {
  try {
    const response = await fetch(`${API}/molecules?limit=100&offset=0`);
    const data = await response.json();
    
    const container = document.getElementById("moleculesTable");
    
    if (data.items.length === 0) {
      container.innerHTML = `
        <div class="empty-state">
          <p>📭 Молекул пока нет. Добавьте первую!</p>
        </div>
      `;
      return;
    }

    container.innerHTML = `
      <p style="color: #666; margin-bottom: 10px;">Всего: ${data.total}</p>
      <table>
        <thead>
          <tr>
            <th>Название</th>
            <th>SMILES</th>
            <th>ID</th>
            <th>Действия</th>
          </tr>
        </thead>
        <tbody>
          ${data.items.map(m => `
            <tr>
              <td>${m.name || "—"}</td>
              <td><span class="smiles-code">${m.smiles}</span></td>
              <td style="font-size: 0.85rem; color: #999;">${m.id.slice(0, 8)}...</td>
              <td>
                <button onclick="deleteMolecule('${m.id}')" class="btn-danger">Удалить</button>
              </td>
            </tr>
          `).join("")}
        </tbody>
      </table>
    `;
  } catch (err) {
    showNotification(`❌ Ошибка загрузки: ${err.message}`, "error");
  }
}

async function deleteMolecule(id) {
  if (!confirm("Удалить эту молекулу?")) return;

  try {
    const response = await fetch(`${API}/molecules/${id}`, { method: "DELETE" });
    if (response.ok) {
      showNotification("✅ Молекула удалена", "success");
      loadMolecules();
    } else {
      showNotification("❌ Ошибка удаления", "error");
    }
  } catch (err) {
    showNotification(`❌ Ошибка: ${err.message}`, "error");
  }
}

async function searchMolecules() {
  const substructure = document.getElementById("substructure").value.trim();
  
  try {
    const response = await fetch(`${API}/search`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ substructure }),
    });

    if (!response.ok) {
      const error = await response.json();
      showNotification(`❌ ${error.detail}`, "error");
      return;
    }

    const data = await response.json();
    const container = document.getElementById("searchResults");

    if (data.matches.length === 0) {
      container.innerHTML = `
        <div class="empty-state">
          <p>🔍 Ничего не найдено</p>
        </div>
      `;
      return;
    }

    container.innerHTML = `
      <p style="margin-top: 20px; color: #667eea; font-weight: 600;">
        Найдено: ${data.matches.length}
      </p>
      <div class="results-grid">
        ${data.matches.map(m => `
          <div class="molecule-card">
            <h3>${m.name || "Без названия"}</h3>
            <div class="smiles">${m.smiles}</div>
            <div class="id">ID: ${m.id.slice(0, 13)}...</div>
          </div>
        `).join("")}
      </div>
    `;
  } catch (err) {
    showNotification(`❌ Ошибка поиска: ${err.message}`, "error");
  }
}

window.addEventListener("DOMContentLoaded", loadMolecules);
